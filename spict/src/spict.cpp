/*
    Stochastic surplus Production model in Continuous-Time (SPiCT)
    Copyright (C) 2015  Martin Waever Pedersen, mawp@dtu.dk or wpsgodd@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// 14.10.2014
#include <TMB.hpp>

/* Predict biomass */
template<class Type>
Type predictlogB(const Type &B0, const Type &F, const Type &gamma, const Type &m, const Type &K, const Type &dt, const Type &n, const Type &sdb2=0)
{
  // Euler discretised Lamperti transformed Pella-Tomlinson surplus production model in Fletcher (1978) form.
  return log(B0) + (gamma*m/K - gamma*m/K*pow(B0/K, n-1.0) - F - 0.5*sdb2)*dt;
}

/* Predict catch*/
template<class Type>
Type predictC(const Type &F, const Type &B0, const Type &dt)
{
  return F*B0*dt;
}

/* Main script */
template<class Type>
Type objective_function<Type>::operator() ()
{
  Type ans=0;

  // DATA
  DATA_INTEGER(reportall);     // Report everything?
  DATA_VECTOR(dt);             // Time steps
  DATA_VECTOR(dtpredcinds);    // Indices of predictions in F state vector 
  DATA_INTEGER(dtpredcnsteps); // Number of sub time step for prediction
  DATA_SCALAR(dtprediind);     // Index of B state vector to use for predicting I
  DATA_INTEGER(indlastobs);    // Index of B and F corresponding to the last observation.
  DATA_VECTOR(obssrt);         // Catch and index observations sorted in time (to enable osar)
  DATA_VECTOR_INDICATOR(keep, obssrt); // This one is required to calculate OSA residuals
  DATA_VECTOR(stdevfacc);      // Factors to scale stdev of catch observation error
  DATA_VECTOR(stdevfaci);      // Factors to scale stdev of index observation error
  DATA_VECTOR(isc);            // Indices in obssrt of catch observations
  DATA_VECTOR(isi);            // Indices in obssrt of index observations
  DATA_INTEGER(nobsC);         // Number of catch observations
  DATA_INTEGER(nobsI);         // Number of index observations
  DATA_VECTOR(ic);             // Vector such that B(ic(i)) is the state at the start of obsC(i)
  DATA_VECTOR(nc);             // nc(i) gives the number of time intervals obsC(i) spans
  DATA_VECTOR(ii);             // A vector such that B(ii(i)) is the state corresponding to I(i)
  DATA_VECTOR(iq);             // A vector such that iq(i) is the index number corresponding to I_iq(i)
  DATA_VECTOR(isdi);           // A vector such that isdi(i) is the index number corresponding to I_isdi(i)
  DATA_VECTOR(ir);             // A vector indicating when the different rs should be used
  DATA_VECTOR(seasons);        // A vector of length ns indicating to which season a state belongs
  DATA_VECTOR(seasonindex);    // A vector of length ns giving the number stepped within the current year
  DATA_MATRIX(splinemat);      // Design matrix for the seasonal spline
  DATA_MATRIX(splinematfine);  // Design matrix for the seasonal spline on a fine time scale to get spline uncertainty
  DATA_SCALAR(omega);          // Period time of seasonal SDEs (2*pi = 1 year period)
  DATA_SCALAR(seasontype);     // Variable indicating whether to use 1=spline, 2=coupled SDEs
  DATA_VECTOR(ffacvec);        // Management factor each year multiply the predicted F with ffac
  DATA_VECTOR(fconvec);        // Management factor each year add this constant to the predicted F
  DATA_VECTOR(indpred);        // A vector indicating when the management factor should be applied
  DATA_SCALAR(robflagc);       // Catch Degrees of freedom of t-distribution (only used if tdf < 25)
  DATA_SCALAR(robflagi);       // Index Degrees of freedom of t-distribution (only used if tdf < 25)
  DATA_INTEGER(stochmsy);      // Use stochastic msy?

  // Priors
  DATA_VECTOR(priorn);         // Prior vector for n, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorr);         // Prior vector for r, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorK);         // Prior vector for K, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorm);         // Prior vector for m, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorq);         // Prior vector for q, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorbkfrac);    // Prior vector for B0/K, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdb);       // Prior vector for sdb, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdf);       // Prior vector for sdf, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdi);       // Prior vector for sdi, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdc);       // Prior vector for sdc, [log(mean), stdev in log, useflag]
  DATA_VECTOR(prioralpha);     // Prior vector for alpha, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorbeta);      // Prior vector for beta, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorB);         // Prior vector for B, [log(mean), stdev in log, useflag, year, ib]
  DATA_VECTOR(priorF);         // Prior vector for F, [log(mean), stdev in log, useflag, year, if]
  DATA_VECTOR(priorBBmsy);     // Prior vector for B/Bmsy, [log(mean), stdev in log, useflag, year, ib]
  DATA_VECTOR(priorFFmsy);     // Prior vector for F/Fmsy, [log(mean), stdev in log, useflag, year, if]

  // Options
  DATA_SCALAR(simple);         // If simple=1 then use simple model (catch assumed known, no F process)
  DATA_SCALAR(dbg);            // Debug flag, if == 1 then print stuff.

  // PARAMETERS
  PARAMETER_VECTOR(logm);      // m following the Fletcher formulation (see Prager 2002)
  PARAMETER(logK);             // Carrying capacity
  PARAMETER_VECTOR(logq);      // Catchability
  PARAMETER(logn);             // Pella-Tomlinson exponent
  PARAMETER(logsdb);           // Standard deviation in B process
  PARAMETER_VECTOR(logsdu);    // Standard deviation in seasonal component of F process
  PARAMETER(logsdf);           // Standard deviation in diffusion component of F process
  PARAMETER_VECTOR(logsdi);    // sdi = alpha*sdb
  PARAMETER(logsdc);           // sdc = beta*sdf
  PARAMETER_VECTOR(logphi);    // Season levels of F.
  PARAMETER(loglambda);        // Damping variable when using seasonal SDEs
  PARAMETER(logitpp);          // Proportion of narrow distribution when using robust obs err.
  PARAMETER(logp1robfac);      // Coefficient to the standard deviation of robust observation error distribution
  PARAMETER_VECTOR(logF);      // Diffusion component of F in log
  PARAMETER_MATRIX(logu);      // Seasonal component of F in log
  PARAMETER_VECTOR(logB);      // Biomass in log


  //std::cout << "expmosc: " << expmosc(lambda, omega, 0.1) << std::endl;
   if(dbg > 0){
     std::cout << "==== DATA read, now calculating derived quantities ====" << std::endl;
   }

  int ind;
  // Distribute sorted observations into logobsC and logobsI vectors
  //int nobsC = isc.size();
  vector<Type> logobsC(nobsC);
  for(int i=0; i<nobsC; i++){ 
    ind = CppAD::Integer(isc(i)-1);
    logobsC(i) = obssrt(ind); 
  }
  //int nobsI = isi.size();
  vector<Type> logobsI(nobsI);
  for(int i=0; i<nobsI; i++){ 
    ind = CppAD::Integer(isi(i)-1);
    logobsI(i) = obssrt(ind);
  }

  vector<Type> logphipar(logphi.size()+1);
  logphipar(0) = 0.0; // The first logphi is set to 0, the rest are estimated relative to this.
  for(int i=1; i<logphipar.size(); i++){ logphipar(i) = logphi(i-1); }
  vector<Type> temp = splinemat.col(0);
  vector<Type> seasonspline(temp.size());
  seasonspline = splinemat * logphipar;
  vector<Type> tempfine = splinematfine.col(0);
  vector<Type> seasonsplinefine(tempfine.size());
  seasonsplinefine = splinematfine * logphipar;
  Type pp = 1.0/(1.0 + exp(-logitpp));
  Type robfac = 1.0 + exp(logp1robfac);
  int nm = logm.size();
  vector<Type> m(nm); 
  for(int i=0; i<nm; i++){ m(i) = exp(logm(i)); }
  Type K = exp(logK);
  int nq = logq.size();
  vector<Type> q(nq);
  for(int i=0; i<nq; i++){ q(i) = exp(logq(i)); }
  vector<Type> logq2(nq);
  for(int i=0; i<nq; i++){ logq2(i) = log(100) + logq(i); }
  Type n = exp(logn);
  Type gamma = pow(n, n/(n-1.0)) / (n-1.0);
  Type p = n - 1.0;
  Type lambda = exp(loglambda);
  Type sdf = exp(logsdf);
  vector<Type> sdu = exp(logsdu);
  int nsdu = sdu.size();
  Type sdb = exp(logsdb);
  Type sdb2 = sdb*sdb;
  int nsdi = logsdi.size();
  //int nalpha = nsdi;
  vector<Type> sdi = exp(logsdi);
  /*
  if(nsdi==1){
    for(int i=0; i<nsdi; i++){ // Same alpha for all indices
      sdi(i) = exp(logsdi(0));
	//alpha(i) = sdi(i)/sdb; 
    } 
  } else {
    sdi = exp(logsdi);
    //alpha = sdi/sdb;
  }
  */
  //vector<Type> alpha(nq);
  vector<Type> alpha = sdi/sdb;
  vector<Type> logalpha = log(alpha);
  int nalpha = logalpha.size();
  Type sdc = exp(logsdc);
  Type beta = sdc/sdf;
  Type logbeta = log(beta);

  // Put wide smooth distributions on difficult parameters to stabilise optimisation
  // Note that this contributes to the objective function, which therefore cannot be regarded as a likelihood
  ans -= dnorm(logB(0), Type(10.0), Type(10.0), 1);
  ans -= dnorm(logF(0), Type(0.0), Type(10.0), 1);
  ans -= dnorm(logn, Type(0.6931472), Type(10.0), 1);
  ans -= dnorm(logbeta, Type(0.0), Type(10.0), 1);
  for(int i=0; i<nalpha; i++){ ans -= dnorm(logalpha(i), Type(0.0), Type(10.0), 1); }


  //vector<Type> logsdi = log(sdi);
  //Type logsdc = log(sdc);
  int nobsCp = ic.size();
  int ns = logF.size();
  vector<Type> P(ns-1);
  vector<Type> B = exp(logB);
  vector<Type> mvec(ns);
  vector<Type> logBmsyvec(ns);
  vector<Type> logFmsyvec(ns);
  //vector<Type> ffacvec(ns);
  //for(int i=0; i<ns; i++){ ffacvec(i) = 1.0; }
  vector<Type> Cpred(nobsCp);
  for(int i=0; i<nobsCp; i++){ Cpred(i) = 0.0; }
  vector<Type> logIpred(nobsI);
  vector<Type> logCpred(nobsCp);

  vector<Type> Bmsyd(nm);
  vector<Type> MSYd = m;
  vector<Type> Fmsyd(nm);
  vector<Type> rbmean(nm);
  vector<Type> Bmsys(nm);
  vector<Type> Fmsys(nm);
  vector<Type> MSYs(nm);
  for(int i=0; i<nm; i++){
    // Deterministic reference points
    Bmsyd(i) = K * pow(1.0/n, 1.0/(n-1.0));
    Fmsyd(i) = MSYd(i)/Bmsyd(i);

    // Stochastic reference points (NOTE: only proved for n>1, Bordet and Rivest (2014))
    rbmean(i) = (n-1)/n*gamma*m(i)/K;
    Bmsys(i) = Bmsyd(i) * (1.0 - (1.0 + rbmean(i)*(p-1.0)/2.0)*sdb2 / (rbmean(i)*pow(2.0-rbmean(i), 2.0)));
    Fmsys(i) = Fmsyd(i) - (p*(1.0-rbmean(i))*sdb2) / pow(2.0-rbmean(i), 2.0) ;
    MSYs(i) = MSYd(i) * (1.0 - ((p+1.0)/2.0*sdb2) / (1.0 - pow(1.0-rbmean(i), 2.0)));
  }

  vector<Type> logBmsyd = log(Bmsyd);
  vector<Type> logMSYd = log(MSYd);
  vector<Type> logFmsyd = log(Fmsyd);
  vector<Type> logBmsys = log(Bmsys);
  vector<Type> logFmsys = log(Fmsys);
  vector<Type> logMSYs = log(MSYs);

  vector<Type> Bmsy(nm);
  vector<Type> MSY(nm);
  vector<Type> Fmsy(nm);
  vector<Type> logBmsy(nm);
  vector<Type> logFmsy(nm);
  vector<Type> logMSY(nm);

  if(stochmsy == 1){
    // Use stochastic reference points
    Bmsy = Bmsys;
    MSY = MSYs;
    Fmsy = Fmsys;
    logBmsy = logBmsys;
    logFmsy = logFmsys;
    logMSY = logMSYs;
  } else {
    // Use deterministic reference points
    Bmsy = Bmsyd;
    MSY = MSYd;
    Fmsy = Fmsyd;
    logBmsy = logBmsyd;
    logFmsy = logFmsyd;
    logMSY = logMSYd;
  }

  vector<Type> Emsy(nq);
  vector<Type> Emsy2(nq);
  for(int i=0; i<nq; i++){ 
    Emsy(i) = Fmsy(0)/q(i); 
    Emsy2(i) = Fmsy(0)/exp(logq2(i))*1.0e4; // Used for the results of the albacore data set
  }
  vector<Type> logEmsy = log(Emsy);
  vector<Type> logEmsy2 = log(Emsy2);

  // Calculate growth rate
  Type sign = 1.0;
  if(n < 1.0) sign = -1.0; // Following Fletcher (1978)
  vector<Type> r(nm);
  vector<Type> logr(nm);
  for(int i=0; i<nm; i++){ 
    //r(i) = sign * gamma * m(i) / K;  // For some reason this doesnt work for n<0
    r(i) =  abs(gamma * m(i) / K);
    logr(i) = log(r(i)); 
    //std::cout << "sign: " << sign << " -- n: " << n << " -- gamma: " << gamma << n << " -- m(i): " << m(i)<< n << " -- K: " << K << " -- r(i): " << r(i) << " -- logr(i): " << logr(i) << std::endl;
  }


  Type likval;

  if(dbg > 0){
    std::cout << "" << std::endl;
    std::cout << "--- DEBUG: script start --- ans: " << ans << std::endl;
    for(int i=0; i<nm; i++){ std::cout << "INPUT: logm(i): " << logm(i) << " -- i: " << i << std::endl; }
    for(int i=0; i<logphi.size(); i++){ std::cout << "INPUT: logphi(i): " << logphi(i) << " -- i: " << i << std::endl; }
    for(int i=0; i<logphipar.size(); i++){ std::cout << "INPUT: logphipar(i): " << logphipar(i) << " -- i: " << i << std::endl; }
    std::cout << "INPUT: logK: " << logK << std::endl;
    for(int i=0; i<nq; i++){ std::cout << "INPUT: logq(i): " << logq(i) << " -- i: " << i << std::endl; }
    std::cout << "INPUT: logn: " << logn << std::endl;
    std::cout << "INPUT: logsdf: " << logsdf << std::endl;
    std::cout << "INPUT: logsdb: " << logsdb << std::endl;
    std::cout << "INPUT: lambda: " << lambda << std::endl;
    std::cout << "INPUT: omega: " << omega << std::endl;
    std::cout << "logobsC.size(): " << logobsC.size() << "  Cpred.size(): " << Cpred.size() << "  logobsI.size(): " << logobsI.size() << "  dt.size(): " << dt.size() << "  logF.size(): " << logF.size() << "  logu.rows(): " << logu.rows() << "  logu.cols(): " << logu.cols() << "  B.size(): " << B.size() << "  P.size(): " << P.size() << "  mvec.size(): " << mvec.size() << "  iq.size(): " << iq.size() << "  ic.size(): " << ic.size() << "  ir.size(): " << ir.size() << "  logFmsy.size(): " << logFmsy.size() << "  logFmsyvec.size(): " << logFmsyvec.size() << "  logBmsy.size(): " << logBmsy.size() << "  logBmsyvec.size(): " << logBmsyvec.size() << "  m.size(): " << m.size() << "  logphi.size(): " << logphi.size() << "  logphipar.size(): " << logphipar.size() << std::endl;
    std::cout << "Bmsys: " << Bmsys << std::endl;
    std::cout << "Fmsys: " << Fmsys << std::endl;
    std::cout << "MSYs: " << MSYs << std::endl;
    std::cout << "Bmsyd: " << Bmsyd << std::endl;
    std::cout << "Fmsyd: " << Fmsyd << std::endl;
    std::cout << "MSYd: " << MSYd << std::endl;
  }
  // Calculate mvec if multiple rs are used (rarely the case).
  for(int i=0; i<ns; i++){
    ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
    if(dbg>1){
      std::cout << "-- i: " << i << " -- ind: " << ind << " -   mvec(i): " << mvec(i) << std::endl;
    }
    mvec(i) = m(ind);
    logFmsyvec(i) = logFmsy(ind);
    logBmsyvec(i) = logBmsy(ind);
  }

  // PRIORS
  if(dbg > 0){
    std::cout << "PRIOR: priorn(0): " << priorn(0) << " -- priorn(1): " << priorn(1) << " -- priorn(2): " << priorn(2) << std::endl;
    std::cout << "PRIOR: priorr(0): " << priorr(0) << " -- priorr(1): " << priorr(1) << " -- priorr(2): " << priorr(2) << std::endl;
    std::cout << "PRIOR: priorK(0): " << priorK(0) << " -- priorK(1): " << priorK(1) << " -- priorK(2): " << priorK(2) << std::endl;
    std::cout << "PRIOR: priorm(0): " << priorm(0) << " -- priorm(1): " << priorm(1) << " -- priorm(2): " << priorm(2) << std::endl;
    std::cout << "PRIOR: priorq(0): " << priorq(0) << " -- priorq(1): " << priorq(1) << " -- priorq(2): " << priorq(2) << std::endl;
  }
  if(priorn(2) == 1) ans-= dnorm(logn, priorn(0), priorn(1), 1); // Prior for logn
  if(priorr(2) == 1 & nm == 1) ans-= dnorm(logr(0), priorr(0), priorr(1), 1); // Prior for logr
  if(priorK(2) == 1) ans-= dnorm(logK, priorK(0), priorK(1), 1); // Prior for logK
  if(priorm(2) == 1 & nm == 1) ans-= dnorm(logm(0), priorm(0), priorm(1), 1); // Prior for logm
  if(priorq(2) == 1 & nq == 1) ans-= dnorm(logq(0), priorq(0), priorq(1), 1); // Prior for logq
  if(priorbkfrac(2) == 1) ans-= dnorm(logB(0) - logK, priorbkfrac(0), priorbkfrac(1), 1); // Prior for logbkfrac
  if(priorsdb(2) == 1) ans-= dnorm(logsdb, priorsdb(0), priorsdb(1), 1); // Prior for logsdb
  if(priorsdf(2) == 1) ans-= dnorm(logsdf, priorsdf(0), priorsdf(1), 1); // Prior for logsdf
  if(priorsdi(2) == 1) for(int i=0; i<nsdi; i++){ ans-= dnorm(logsdi(i), priorsdi(0), priorsdi(1), 1); } // Prior for logsdi
  if(priorsdc(2) == 1) ans-= dnorm(logsdc, priorsdc(0), priorsdc(1), 1); // Prior for logsdc
  if(prioralpha(2) == 1) for(int i=0; i<nalpha; i++){ ans-= dnorm(logalpha(i), prioralpha(0), prioralpha(1), 1); } // Prior for logalpha
  if(priorbeta(2) == 1) ans-= dnorm(logbeta, priorbeta(0), priorbeta(1), 1); // Prior for logbeta
  if(priorB(2) == 1){
    ind = CppAD::Integer(priorB(4)-1);
    ans-= dnorm(logB(ind), priorB(0), priorB(1), 1); // Prior for logB
  }
  if(priorF(2) == 1){
    ind = CppAD::Integer(priorF(4)-1);
    ans-= dnorm(logF(ind), priorF(0), priorF(1), 1); // Prior for logF
  }
  if(priorBBmsy(2) == 1){
    ind = CppAD::Integer(priorBBmsy(4)-1);
    ans-= dnorm(logB(ind) - logBmsyvec(ind), priorBBmsy(0), priorBBmsy(1), 1); // Prior for logBBmsy
  }
  if(priorFFmsy(2) == 1){
    ind = CppAD::Integer(priorFFmsy(4)-1);
    ans-= dnorm(logF(ind) - logFmsyvec(ind), priorFFmsy(0), priorFFmsy(1), 1); // Prior for logFFmsy
  }

  //for(int i=1; i<indpred.size(); i++){ // don't use i=0 because this is only for plotting
  //  ind = CppAD::Integer(indpred(i)-1); // minus 1 because R starts at 1 and c++ at 0
  //  if(dbg>1){
  //    std::cout << "-- i: " << i << " -- ind: " << ind << " -   ffacvec(i): " << ffacvec(i) << std::endl;
  //  }
  //  ffacvec(ind) = ffac;
  //}

  /*
  dt[i] is the length of the time interval between t_i and t_i+1
  B[t] is biomass at the beginning of the time interval starting at time t
  Binf[t] is equilibrium biomass with the parameters given at the beginning of the time interval starting at time t
  I[t] is an index of biomass (e.g. CPUE) at the beginning of the time interval starting at time t
  P[t] is the accumulated biomass production over the interval starting at time t
  F[t] is the constant fishing mortality during the interval starting at time t
  rvec[t] is the constant intrinsic growth rate during the interval starting at time t
  C[t] is the catch removed during the interval starting at time t.
  ffacvec[t] is the factor to multiply F[t-1] when calculating F[t] such that management is imposed at time t.
  */

  /*
  --- PROCESS EQUATIONS ---
  */

  // FISHING MORTALITY
  vector<Type> logFs = logF;
  if(simple==0){
    if(dbg>0){
      std::cout << "--- DEBUG: F loop start --- ans: " << ans << std::endl;
    }
    // Diffusion component of F
    for(int i=1; i<ns; i++){
      Type logFpred = log( ffacvec(i) * exp(logF(i-1)) + fconvec(i) );
      likval = dnorm(logF(i), logFpred, sqrt(dt(i-1))*sdf, 1);
      ans-=likval;
      // DEBUGGING
      if(dbg>1){
	std::cout << "-- i: " << i << " -   logF(i-1): " << logF(i-1) << "  logF(i): " << logF(i) << "  ffacvec(i): " << ffacvec(i) << "  fconvec(i): " << fconvec(i) << "  sdf: " << sdf << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }

    // Seasonal components
    if(dbg>0){ std::cout << "-- seasontype: " << seasontype << std::endl; }
    //vector<Type> logFs(ns);
    if(seasontype == 1.0){
      // Spline
      int ind2;
      for(int i=0; i<ns; i++){
	ind2 = CppAD::Integer(seasonindex(i));
	logFs(i) += seasonspline(ind2);
	// DEBUGGING
	if(dbg>1){
	  std::cout << "-- i: " << i << " -   logF(i): " << logF(i) << " logFs(i): " << logFs(i) << " ind2: " << ind2 << " seasonspline(ind2): " << seasonspline(ind2) << std::endl;
	}
      }
    }
    if(seasontype == 2.0){
      // Coupled SDEs
      for(int j=0; j<nsdu; j++){
	Type per = j+1.0;
	if(dbg>0){ std::cout << "-- j:" << j << "- per:" << per << "- omega:" << omega << std::endl; } 
	for(int i=1; i<ns; i++){
	  // Analytical expression for matrix exponential
	  matrix<Type> trigmat(2, 2);
	  trigmat(0, 0) = cos(per*omega*dt(i-1));
	  trigmat(0, 1) = -sin(per*omega*dt(i-1));
	  trigmat(1, 0) = sin(per*omega*dt(i-1));
	  trigmat(1, 1) = cos(per*omega*dt(i-1));
	  if(dbg>0){ std::cout << "-- trigmat: " << trigmat << std::endl; }
	  matrix<Type> expmAt = exp(-lambda*dt(i-1))*trigmat;
	  if(dbg>0){ std::cout << "-- expmAt: " << expmAt << std::endl; }
	  // Corrected standard deviation
	  Type sduana = sdu(j) * sqrt(1.0/(2.0*lambda) * (1.0 - exp(-2.0*lambda*dt(i-1))));
	  if(dbg>0){ std::cout << "-- sduana: " << sduana << std::endl; }
	  vector<Type> sublogumF = logu.col(i-1);
	  vector<Type> sublogum(2);
	  sublogum(0) = sublogumF(2*j);
	  sublogum(1) = sublogumF(2*j+1);
	  if(dbg>0){ std::cout << "-- sublogumF: " << sublogumF << "-- sublogum: " << sublogum << std::endl; }
	  vector<Type> logupred = expmAt * sublogum;
	  if(dbg>0){ std::cout << "-- logupred: " << logupred << std::endl; }
	  likval = 0.0;
	  for(int k=0; k<logupred.size(); k++){ 
	    if(dbg>0){ std::cout << "-- k: " << k << "- 2*j+k: " << 2*j+k << " - logu(2*j+k, i): " << logu(2*j+k, i) << std::endl; }
	    likval += dnorm(logu(2*j+k, i), logupred(k), sduana, 1); 
	  }
	  ans-=likval;
	  // DEBUGGING
	  if(dbg>0){
	    std::cout << "-- i: " << i << " -   logu(0,i): " << logu(0,i) << "  logupred(0): " << logupred(0) << " -   logu(1,i): " << logu(1,i) << "  logupred(1): " << logupred(1) << "  sdu(j): " << sdu(j) << "  sduana: " << sduana << "  likval: " << likval << "  ans:" << ans << std::endl;
	  }
	}
	for(int i=0; i<ns; i++) logFs(i) += logu(2*j, i); // Sum diffusion and seasonal component
      }
    }
  } else {
    for(int i=0; i<ns; i++) logFs(i) = -30; // If using simple set fishing mortality to something small.
  }
  vector<Type> F = exp(logFs);

  // BIOMASS PREDICTIONS
  if(dbg>0){
    std::cout << "--- DEBUG: B loop start --- ans: " << ans << std::endl;
  }
  vector<Type> logBpred(ns);
  for(int i=0; i<(ns-1); i++){
    // To predict B(i) use dt(i-1), which is the time interval from t_i-1 to t_i
    if(simple==0){
      logBpred(i+1) = predictlogB(B(i), F(i), gamma, mvec(i), K, dt(i), n, sdb2);
    } else {
      Type Ftmp = 0.0;
      Type Bpredtmp = exp(predictlogB(B(i), Ftmp, gamma, mvec(i), K, dt(i), n, sdb2)) - exp(logobsC(i));
      if(Bpredtmp < 0) Bpredtmp = 1e-8; // Ugly ugly ugly hack to avoid taking log of negative
      logBpred(i+1) = log(Bpredtmp);
      logFs(i) = logobsC(i) - logB(i); // Calculate fishing mortality
    }
    likval = dnorm(logBpred(i+1), logB(i+1), sqrt(dt(i))*sdb, 1);
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logB(i+1): " << logB(i+1) << "  logBpred(i+1): " << logBpred(i+1) << "  sdb: " << sdb << "  likval: " << likval << "  ans:" << ans << std::endl;
    }
  }
  if(simple==1){ logFs(ns-1) = logFs(ns-2);}

  // CATCH PREDICTIONS
  vector<Type> Cpredsub(ns);
  if(simple==0){
    for(int i=0; i<(ns-1); i++){ // ns-1 because dt is 1 shorter than state vec
      // For Cpredsub(i) use dt(i) because Cpredsub(i) is integrated over t_i to t_i+1
      Cpredsub(i) =  predictC(F(i), B(i), dt(i));
    }
    for(int i=0; i<nobsCp; i++){
      // Sum catch contributions from each sub interval
      for(int j=0; j<nc(i); j++){
	ind = CppAD::Integer(ic(i)-1) + j; // minus 1 because R starts at 1 and c++ at 0
	Cpred(i) += Cpredsub(ind);
      }
      logCpred(i) = log(Cpred(i));
      // DEBUGGING
      if(dbg>1){
	std::cout << "-- i: " << i << " -  ind: " << ind << "  logCpred(i): " << logCpred(i) << std::endl;
      }
    }
  } else {
    for(int i=0; i<(ns-1); i++){ // ns-1 because logobsC is 1 shorter than ns
      Cpredsub(i) = exp(logobsC(i));
      logCpred(i) = logobsC(i);
    }
  }

  // CALCULATE PRODUCTION
  for(int i=0; i<(ns-1); i++) P(i) = B(i+1) - B(i) + Cpredsub(i);


  /*
    --- OBSERVATION EQUATIONS ---
  */
  int inds;
  // CATCHES
  if(simple==0){
    if(dbg>0){
      std::cout << "--- DEBUG: Cpred loop start --- ans: " << ans << std::endl;
    }
    // fac and pp are used for the outlier robust Gaussian mixture.
    for(int i=0; i<nobsC; i++){
      //int j = CppAD::Integer(nc(i)-1); <-- NOT USED?
      //ind = CppAD::Integer(ic(i)-1) + j; // minus 1 because R starts at 1 and c++ at 0 <- NOT USED?
      inds = CppAD::Integer(isc(i)-1);
      if(robflagc==1.0){
	likval = log(pp*dnorm(logCpred(i), logobsC(i), stdevfacc(i)*sdc, 0) + (1.0-pp)*dnorm(logCpred(i), logobsC(i), robfac*stdevfacc(i)*sdc, 0));
      } else {
	likval = dnorm(logCpred(i), logobsC(i), stdevfacc(i)*sdc, 1);
      }
      ans-= keep(inds) * likval;
      // DEBUGGING
      if(dbg>1){
	std::cout << "-- i: " << i << " -   logobsC(i): " << logobsC(i) << " -   stdevfacc(i): " << stdevfacc(i) << "  sdc: " << sdc << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
  }

  // BIOMASS INDEX
  if(dbg>0){
    std::cout << "--- DEBUG: Ipred loop start --- ans: " << ans << std::endl;
    std::cout << " logobsI.size(): " << logobsI.size() << "  iq.size(): " << iq.size() << "  ii.size(): " << ii.size() << "  logq.size(): " << logq.size() << "  nobsI: " << nobsI << std::endl;
  }
  int indq;
  int indsdi;
  for(int i=0; i<nobsI; i++){
    ind = CppAD::Integer(ii(i)-1);
    indq = CppAD::Integer(iq(i)-1);
    indsdi = CppAD::Integer(isdi(i)-1);
    inds = CppAD::Integer(isi(i)-1);
    logIpred(i) = logq(indq) + log(B(ind));
    if(robflagi==1.0){
      likval = log(pp*dnorm(logobsI(i), logIpred(i), stdevfaci(i)*sdi(indsdi), 0) + (1.0-pp)*dnorm(logobsI(i), logIpred(i), robfac*stdevfaci(i)*sdi(indsdi), 0));
    } else {
      likval = dnorm(logobsI(i), logIpred(i), stdevfaci(i)*sdi(indsdi), 1);
    }
    ans-= keep(inds) * likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -  ind: " << ind << " -  indq: " << indq << " -  indsdi: " << indsdi << " -  inds: " << inds << " -   logobsI(i): " << logobsI(i) << "  logIpred(i): " << logIpred(i) << "  stdevfaci(i): " << stdevfaci(i) << "  likval: " << likval << "  sdi: " << sdi << "  ans:" << ans << std::endl;
    }
  }

  /*
  --- ONE-STEP-AHEAD PREDICTIONS ---
  */

  if(dbg>0){
    std::cout << "--- DEBUG: ONE-STEP-AHEAD PREDICTIONS --- ans: " << ans << std::endl;
    std::cout << "-- dtpredcnsteps: " << dtpredcnsteps << "  dtpredcinds.size(): " << dtpredcinds.size() <<std::endl;
  }
  Type Cp = 0.0;
  for(int i=0; i<dtpredcnsteps; i++){
    ind = CppAD::Integer(dtpredcinds(i)-1);
    if(dbg>1){
      std::cout << "-- i: " << i << " -  dtpredcinds(i)-1: " << ind << std::endl;
    }
    Cp += Cpredsub(ind);
  }
  Type logCp = log(Cp);

  // Biomass and F at the end of the prediction time interval
  int pind = CppAD::Integer(dtprediind-1);
  Type Bp = B(pind); 
  Type logBp = log(Bp);
  Type logBpBmsy = logBp - logBmsyvec(pind);
  Type logBpK = logBp - logK;
  Type logFp = logFs(pind); 
  Type logFpFmsy = logFp - logFmsyvec(pind);

  vector<Type> logIp(nq);
  for(int i=0; i<nq; i++){
    logIp(i) = logq(i) + log(Bp);
  }

  // Biomass and fishing mortality at last time point
  Type logBl = logB(indlastobs-1);
  Type logBlBmsy = logBl - logBmsyvec(indlastobs-1);
  Type logBlK = logBl - logK;
  Type logFl = logFs(indlastobs-1);
  Type logFlFmsy = logFl - logFmsyvec(indlastobs-1);

  // Calculate relative levels of biomass and fishing mortality
  vector<Type> logBBmsy(ns);
  vector<Type> logFFmsy(ns);
  for(int i=0; i<ns; i++){ 
    logBBmsy(i) = logB(i) - logBmsyvec(i); 
    logFFmsy(i) = logFs(i) - logFmsyvec(i); 
  }

  // 
  Type logbkfrac = logB(0) - logK;

  // ADREPORTS
  ADREPORT(Bmsy);  
  ADREPORT(Bmsyd);
  ADREPORT(Bmsys);
  ADREPORT(logBmsy);
  ADREPORT(logBmsyd);
  ADREPORT(logBmsys);
  ADREPORT(logBp);
  ADREPORT(logBpBmsy);
  ADREPORT(logBpK);
  ADREPORT(logBl);
  ADREPORT(logBlBmsy);
  ADREPORT(logBlK);
  ADREPORT(Fmsy);
  ADREPORT(Fmsyd);
  ADREPORT(Fmsys);
  ADREPORT(logFmsy);
  ADREPORT(logFmsyd);
  ADREPORT(logFmsys);
  ADREPORT(logFp);
  ADREPORT(logFpFmsy);
  ADREPORT(logFl);
  ADREPORT(logFlFmsy);
  ADREPORT(MSY);
  ADREPORT(MSYd);
  ADREPORT(MSYs);
  ADREPORT(logMSY);
  ADREPORT(logMSYd);
  ADREPORT(logMSYs);
  ADREPORT(Emsy);
  ADREPORT(Emsy2);
  ADREPORT(logEmsy);
  ADREPORT(logEmsy2);
  ADREPORT(logbkfrac);
  ADREPORT(seasonsplinefine);
  // PREDICTIONS
  ADREPORT(Cp);
  ADREPORT(logIp);
  ADREPORT(logCp);
  // PARAMETERS
  ADREPORT(r);
  ADREPORT(logr);
  ADREPORT(K);
  ADREPORT(q);
  //ADREPORT(logq);
  ADREPORT(logq2);
  ADREPORT(p);
  ADREPORT(gamma);
  ADREPORT(m);
  ADREPORT(sdf);
  ADREPORT(sdc);
  ADREPORT(sdb);
  ADREPORT(sdi);
  ADREPORT(logalpha);
  ADREPORT(logbeta);
  if(reportall){ 
    // These reports are derived from the random effects and are therefore vectors. TMB calculates the covariance of all sdreports leading to a very large covariance matrix which may cause memory problems.
    // B
    ADREPORT(logBBmsy);
    // F
    ADREPORT(logFFmsy); // Vector of size ns
    ADREPORT(logFs);    // Vector of size ns
    // C
    ADREPORT(logCpred);
    // Other
    ADREPORT(logIpred);
  }

  // REPORTS (these don't require sdreport to be output)
  REPORT(Cp);
  REPORT(logIp);
  REPORT(MSY);
  REPORT(Bmsy);
  REPORT(Fmsy);
  REPORT(stochmsy);

  return ans;
}

