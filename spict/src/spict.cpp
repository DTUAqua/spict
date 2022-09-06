/*
  Stochastic surplus Production model in Continuous-Time (SPiCT)
  Copyright (C) 2015-2016  Martin W. Pedersen, mawp@dtu.dk, wpsgodd@gmail.com

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
#define TMB_LIB_INIT R_init_spict
#include <TMB.hpp>

/* Predict biomass */
template<class Type>
Type predictlogB(const Type &B0, const Type &F, const Type &gamma, const Type &m, const Type &K, const Type &dt, const Type &n, const Type &sdb2)
{
  // Euler discretised Lamperti transformed Pella-Tomlinson surplus production model in Fletcher (1978) form.
  return log(B0) + (gamma*m/K - gamma*m/K*pow(B0/K, n-1.0) - F - 0.5*sdb2)*dt;
}

/* Predict m */
template<class Type>
Type predictm(const Type &logm0, const Type &dt, const Type &sdm2, const Type &psi)
{
  //return logm0 + psi*(logmc - logm0)*dt;
  return logm0 - psi*logm0*dt;
}

/* Predict F1 */
template<class Type>
Type predictF1(const Type &logF0, const Type &dt, const Type &sdf2, const Type &delta, const Type &logeta)
{
  return exp(logF0 + delta*(logeta - logF0)*dt);
}

/* Predict F2 */
template<class Type>
Type predictF2(const Type &logF0, const Type &dt, const Type &sdf2, const Type &delta, const Type &logeta)
{
  return exp(logF0 - 0.5*sdf2*dt + delta*(logeta - logF0)*dt);
}

/* Predict catch */
template<class Type>
Type predictC(const Type &F, const Type &B0, const Type &dt)
{
  return F*B0*dt;
}

/* Step function, 0 if n < 1, 1 if n > 1 */
// Note that p = n - 1;
template<class Type>
Type stepfun(const Type &p)
{
  return 1.0 / (1.0 + exp(-10000 * p));
}

template <class Type>
Type ilogit(Type x){
  return Type(1.0)/(Type(1.0)+exp(-x));
}


/* Main script */
template<class Type>
Type objective_function<Type>::operator() ()
{
  Type ans=0;

  // DATA
  DATA_INTEGER(reportall);     // Report everything?
  DATA_INTEGER(reportRel);     // ADreport B/mean(B) and F/mean(F) ?
  DATA_VECTOR(dt);             // Time steps
  DATA_VECTOR(dtpredcinds);    // Indices of predictions in F state vector
  DATA_INTEGER(dtpredcnsteps); // Number of sub time step for prediction
  DATA_SCALAR(dtprediind);     // Index of B state vector to use for predicting I
  DATA_VECTOR(dtpredeinds);    // Indices of predictions in F state vector
  DATA_INTEGER(dtpredensteps); // Number of sub time step for prediction
  DATA_INTEGER(indlastobs);    // Index of B and F corresponding to the last observation.
  DATA_VECTOR(obssrt);         // Catch and index observations sorted in time (to enable osar)
  DATA_VECTOR_INDICATOR(keep, obssrt); // This one is required to calculate OSA residuals
  DATA_VECTOR(stdevfacc);      // Factors to scale stdev of catch observation error
  DATA_VECTOR(stdevfaci);      // Factors to scale stdev of index observation error
  DATA_VECTOR(stdevface);      // Factors to scale stdev of effort observation error
  DATA_VECTOR(isc);            // Indices in obssrt of catch observations
  DATA_VECTOR(isi);            // Indices in obssrt of index observations
  DATA_VECTOR(ise);            // Indices in obssrt of effort observations
  DATA_INTEGER(nobsC);         // Number of catch observations
  DATA_INTEGER(nobsI);         // Number of index observations
  DATA_INTEGER(nobsE);         // Number of effort observations
  DATA_VECTOR(ic);             // Vector such that B(ic(i)) is the state at the start of obsC(i)
  DATA_VECTOR(nc);             // nc(i) gives the number of time intervals obsC(i) spans
  DATA_VECTOR(ie);             // Vector such that E(ie(i)) is the state at the start of obsE(i)
  DATA_VECTOR(ne);             // ne(i) gives the number of time intervals obsE(i) spans
  DATA_VECTOR(ii);             // A vector such that B(ii(i)) is the state corresponding to I(i)
  DATA_VECTOR(iq);             // A vector such that iq(i) is the index number corresponding to I_iq(i)
  DATA_VECTOR(isdi);           // A vector such that isdi(i) is the index number corresponding to I_isdi(i)
  DATA_VECTOR(ir);             // A vector indicating when the different rs should be used
  DATA_VECTOR(isdf);           //
  DATA_VECTOR(logmcov);        // A vector containing covariate information for logm
  DATA_VECTOR(seasons);        // A vector of length ns indicating to which season a state belongs
  DATA_VECTOR(seasonindex);    // A vector of length ns giving the number stepped within the current year
  DATA_INTEGER(nseasons);      // Number of seasons pr year
  DATA_VECTOR(seasonindex2)    // A vector of length ns mapping states to seasonal AR component (for seasontype=3)
    DATA_MATRIX(splinemat);      // Design matrix for the seasonal spline
  DATA_MATRIX(splinematfine);  // Design matrix for the seasonal spline on a fine time scale to get spline uncertainty
  DATA_SCALAR(omega);          // Period time of seasonal SDEs (2*pi = 1 year period)
  DATA_INTEGER(seasontype);     // Variable indicating whether to use 1=spline, 2=coupled SDEs
  DATA_INTEGER(efforttype);     // Variable indicating whether to use 1=spline, 2=coupled SDEs
  DATA_INTEGER(timevaryinggrowth); //  Flag indicating whether REs are used for growth
  DATA_INTEGER(logmcovflag);   // Flag indicating whether covariate information is available
  DATA_VECTOR(ffacvec);        // Management factor each year multiply the predicted F with ffac
  DATA_VECTOR(fconvec);        // Management factor each year add this constant to the predicted F
  DATA_INTEGER(robflagc);       // If 1 use robust observation error for catches
  DATA_IVECTOR(robflagi);       // If 1 use robust observation error for index
  DATA_INTEGER(robflage);       // If 1 use robust observation error for effort
  DATA_INTEGER(stochmsy);      // Use stochastic msy?
  DATA_INTEGER(stabilise);     // If 1 stabilise optimisation using uninformative priors
  //DATA_SCALAR(effortflag);     // If effortflag == 1 use effort data, else use index data
  DATA_FACTOR(MSYregime);      // factor mapping each time step to an m-regime
  DATA_VECTOR(iuse);

  // Priors
  DATA_VECTOR(priorn);         // Prior vector for n, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorngamma);    // Prior vector for logn, gamma distribution [shape, rate, useflag ]
  DATA_VECTOR(priorr);         // Prior vector for r, [log(mean), stdev in log, useflag]
  //DATA_VECTOR(priorrp);        // Prior vector for rp, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorK);         // Prior vector for K, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorm);         // Prior vector for m, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priormu);        // Prior vector for mu, [log(mean), stdev in log, useflag]
  DATA_MATRIX(priorq);         // Prior vector for q, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorqf);        // Prior vector for qf, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorbkfrac);    // Prior vector for B0/K, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdb);       // Prior vector for sdb, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdm);       // Prior vector for sdm, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdf);       // Prior vector for sdf, [log(mean), stdev in log, useflag]
  DATA_MATRIX(priorsdi);       // Prior vector for sdi, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsde);       // Prior vector for sde, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdc);       // Prior vector for sdc, [log(mean), stdev in log, useflag]
  DATA_VECTOR(prioralpha);     // Prior vector for alpha, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorbeta);      // Prior vector for beta, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorpsi);       // Prior vector for psi, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorB);         // Prior vector for B, [log(mean), stdev in log, useflag, year, ib]
  DATA_VECTOR(priorF);         // Prior vector for F, [log(mean), stdev in log, useflag, year, if]
  DATA_VECTOR(priorBBmsy);     // Prior vector for B/Bmsy, [log(mean), stdev in log, useflag, year, ib]
  DATA_VECTOR(priorFFmsy);     // Prior vector for F/Fmsy, [log(mean), stdev in log, useflag, year, if]
  DATA_VECTOR(priorBmsyB0)     // Prior vector for Bmsy/B_0, [mean, stdev, useflag]
    // Options
    DATA_SCALAR(simple);         // If simple=1 then use simple model (catch assumed known, no F process)
  DATA_SCALAR(dbg);            // Debug flag, if == 1 then print stuff.
  DATA_INTEGER(reportmode);    // If 1-5 only specific quantities are ADreported (increases speed, relevant for fitting within MSE)
  DATA_INTEGER(simRandomEffects); // flag turning simulation of random effects on/off

  // PARAMETERS
  PARAMETER_VECTOR(logm);      // m following the Fletcher formulation (see Prager 2002)
  PARAMETER(mu);               // Coefficient for covariate info for logm
  PARAMETER(logK);             // Carrying capacity
  PARAMETER_VECTOR(logq);      // Catchability for index
  PARAMETER(logqf);            // Catchability for effort
  PARAMETER(logn);             // Pella-Tomlinson exponent
  PARAMETER(logsdb);           // Standard deviation in B process
  PARAMETER_VECTOR(logsdu);    // Standard deviation in seasonal component of F process
  PARAMETER_VECTOR(logsdf);    // Standard deviation in diffusion component of F process
  PARAMETER_VECTOR(logsdi);    // sdi = alpha*sdb
  PARAMETER(logsde);           // sdc = beta*sdf
  PARAMETER(logsdc);           // sdc = beta*sdf
  PARAMETER(logsdm);           //
  PARAMETER(logpsi);           // Mean reversion in OU for logm
  PARAMETER_VECTOR(logphi);    // Season levels of F.
  PARAMETER(loglambda);        // Damping variable when using seasonal SDEs
  PARAMETER(logdelta);          // Strength of mean reversion in OU F process (delta = 0 mean RW)
  PARAMETER(logeta);           // Mean of OU F process
  PARAMETER(logitpp);          // Proportion of narrow distribution when using robust obs err.
  PARAMETER(logp1robfac);      // Coefficient to the standard deviation of robust observation error distribution
  PARAMETER_VECTOR(logF);      // Diffusion component of F in log
  PARAMETER_MATRIX(logu);      // Seasonal component of F in log
  PARAMETER_VECTOR(logB);      // Biomass in log
  PARAMETER_VECTOR(logmre);    // Random effect on m
  PARAMETER_VECTOR(SARvec);    // Autoregressive deviations to seasonal spline
  PARAMETER(logitSARphi);      // AR coefficient for seasonal spline dev
  PARAMETER(logSdSAR);         // Standard deviation seasonal spline deviations


  //std::cout << "expmosc: " << expmosc(lambda, omega, 0.1) << std::endl;
  if(dbg > 0){
    std::cout << "==== DATA read, now calculating derived quantities ====" << std::endl;
  }

  int ind = 0;
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
  vector<Type> logobsE(nobsE);
  for(int i=0; i<nobsE; i++){
    ind = CppAD::Integer(ise(i)-1);
    logobsE(i) = obssrt(ind);
  }

  // Length of vectors
  int nm = logm.size();
  int nq = logq.size();
  int nsdf = logsdf.size();
  int nsdu = logsdu.size();
  int nsdi = logsdi.size();
  int nobsCp = ic.size();
  int ns = logF.size();

  // Transform parameters
  Type psi = exp(logpsi);
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
  vector<Type> m(nm);
  for(int i=0; i<nm; i++){ m(i) = exp(logm(i)); }
  Type K = exp(logK);
  vector<Type> q(nq);
  for(int i=0; i<nq; i++){ q(i) = exp(logq(i)); }
  vector<Type> logq2(nq);
  for(int i=0; i<nq; i++){ logq2(i) = log(100) + logq(i); }
  //Type qf = exp(logqf);

  Type lambda = exp(loglambda);
  Type delta = exp(logdelta);
  vector<Type> sdf(nsdf);
  vector<Type> sdf2(nsdf);
  vector<Type> isdf2(nsdf);
  for(int i=0; i<nsdf; i++){
    sdf(i) = exp(logsdf(i));
    sdf2(i) = sdf(i)*sdf(i);
    isdf2(i) = 1.0/sdf2(i);

  }
  vector<Type> sdu = exp(logsdu);
  Type sdb = exp(logsdb);
  Type sdb2 = sdb*sdb;
  Type isdb2 = 1.0/sdb2;
  Type sdm = exp(logsdm);
  Type sdm2 = sdm*sdm;
  //Type isdm2 = 1.0/sdm2;
  Type sde = exp(logsde);
  Type sde2 = sde*sde;
  Type isde2 = 1.0/sde2;
  Type SARphi = ilogit(logitSARphi);
  Type sdSAR = exp(logSdSAR);


  // Initialise vectors
  vector<Type> P(ns-1);
  vector<Type> B = exp(logB);
  //vector<Type> mvec(ns);
  vector<Type> Cpred(nobsCp);
  for(int i=0; i < nobsCp; i++){
    Cpred(i) = 0.0;
  }
  vector<Type> logIpred(nobsI);
  vector<Type> logCpred(nobsCp);
  vector<Type> logEpred(nobsE);

  // Covariate for m
  vector<Type> logmc(ns);
  for(int i=0; i < ns; i++){
    logmc(i) = logm(MSYregime[i]) + mu*logmcov(i);
  }

  // Reference points
  vector<Type> mvec(ns);
  for(int i=0; i < ns; i++){
    //mvec(i) = exp(logm(0) + mu*logmcov(i) + logmre(i));
    mvec(i) = exp(logmc(i) + logmre(i));
  }

  // Parameters with default priors
  Type n = exp(logn);
  Type gamma = pow(n, n/(n-1.0)) / (n-1.0);
  vector<Type> sdi = exp(logsdi);
  vector<Type> sdi2(nsdi);
  vector<Type> isdi2(nsdi);
  vector<Type> alpha = sdi/sdb;
  vector<Type> logalpha = log(alpha);
  for(int i=0; i<nsdi; i++){
    sdi2(i) = sdi(i)*sdi(i);
    isdi2(i) = 1.0/sdi2(i);
  }
  Type sdc = exp(logsdc);
  Type sdc2 = sdc*sdc;
  Type isdc2 = 1.0/sdc2;
  Type beta = sdc/sdf(0);
  Type logbeta = log(beta);


  Type p = n - 1.0;
  vector<Type> Bmsyd(nm);
  vector<Type> Fmsyd(nm);
  vector<Type> MSYd = m;
  vector<Type> Bmsys(nm);
  vector<Type> Fmsys(nm);
  vector<Type> MSYs(nm);
  int flag = asDouble(n) > 1; // Cast n as double to calc flag
  for(int i=0; i<nm; i++){
    // Deterministic reference points
    Bmsyd(i) = K * pow(1.0/n, 1.0/(n-1.0));
    Fmsyd(i) = MSYd(i)/Bmsyd(i);
    // Stochastic reference points (NOTE: only proved for n>1, Bordet and Rivest (2014))
    // The stepfun ensures that stochastic reference points are only used if n > 1.
    Type BmsyStochContr = (1.0 - (1.0 + Fmsyd(i)*(p-1.0)/2.0)*sdb2 / (Fmsyd(i)*pow(2.0-Fmsyd(i), 2.0)));
    Type FmsyStochContr = (p*(1.0-Fmsyd(i))*sdb2) / pow(2.0-Fmsyd(i), 2.0);
    Type MSYstochContr = (1.0 - ((p+1.0)/2.0*sdb2) / (1.0 - pow(1.0-Fmsyd(i), 2.0)));

    //flag = asDouble(n) > 1 & asDouble(BmsyStochContr) > 0;
    if(flag){
      Bmsys(i) = Bmsyd(i) * BmsyStochContr;
      Fmsys(i) = Fmsyd(i) - FmsyStochContr;
      MSYs(i) = MSYd(i) * MSYstochContr;
    } else {
      Bmsys(i) = Bmsyd(i);
      Fmsys(i) = Fmsyd(i);
      MSYs(i) = MSYd(i);
    }
    //std::cout << "flag: " << flag << std::endl;
    //std::cout << "BmsyStochContr: " << BmsyStochContr << std::endl;
    //std::cout << "Bmsyd(i): " << Bmsyd(i) << std::endl;
    //std::cout << "Bmsys(i): " << Bmsys(i) << std::endl;
  }

  // log reference points
  vector<Type> logBmsyd = log(Bmsyd);
  vector<Type> logMSYd = log(MSYd);
  vector<Type> logFmsyd = log(Fmsyd);
  vector<Type> logBmsys = log(Bmsys);
  vector<Type> logFmsys = log(Fmsys);
  vector<Type> logMSYs = log(MSYs);
  // Used reference points
  vector<Type> Bmsy(nm);
  vector<Type> MSY(nm);
  vector<Type> Fmsy(nm);
  vector<Type> logBmsy(nm);
  vector<Type> logFmsy(nm);
  vector<Type> logMSY(nm);
  // Reference point vectors (when time varying growth)
  vector<Type> logFmsyvec(ns);
  vector<Type> logBmsyvec(ns);
  vector<Type> logMSYvec(ns);

  vector<Type> Bmsy2(nm);
  if(flag){
    Bmsy2 = Bmsys;
  } else {
    Bmsy2 = Bmsyd;
  }

  if(stochmsy == 1){
    // Use stochastic reference points
    Bmsy = Bmsys;
    MSY = MSYs;
    Fmsy = Fmsys;
    logBmsy = logBmsys;
    logFmsy = logFmsys;
    logMSY = logMSYs;
    for(int i=0; i < ns; i++){
      ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
      Type Fmsydveci = mvec(i) / Bmsyd(ind);
      logFmsyvec(i) = log(Fmsydveci - (p*(1.0-Fmsydveci)*sdb2) / pow(2.0-Fmsydveci, 2.0));
      logBmsyvec(i) = logBmsys(ind);
      logMSYvec(i) = log(mvec(i) * (1.0 - ((p+1.0)/2.0*sdb2) / (1.0 - pow(1.0-Fmsydveci, 2.0))));
    }
  } else {
    // Use deterministic reference points
    Bmsy = Bmsyd;
    MSY = MSYd;
    Fmsy = Fmsyd;
    logBmsy = logBmsyd;
    logFmsy = logFmsyd;
    logMSY = logMSYd;
    for(int i=0; i<ns; i++){
      ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
      logFmsyvec(i) = log(mvec(i) / Bmsyd(ind));
      logBmsyvec(i) = logBmsyd(ind);
      logMSYvec(i) = log(mvec(i));
    }
  }

  // These quantities are calculated to enable comparison with the Polacheck et al (1993) parameter estimates
  vector<Type> Emsy(nq);
  vector<Type> Emsy2(nq);
  for(int i=0; i<nq; i++){
    Emsy(i) = Fmsy(0)/q(i);
    Emsy2(i) = Fmsy(0)/exp(logq2(i))*1.0e4; // Used for the results of the albacore data set
  }
  vector<Type> logEmsy = log(Emsy);
  vector<Type> logEmsy2 = log(Emsy2);

  // Calculate growth rate
  //Type sign = 1.0;
  //if(n < 1.0) sign = -1.0; // Following Fletcher (1978)
  vector<Type> r(nm);
  vector<Type> logr(nm);
  vector<Type> logrre(ns);
  vector<Type> rc(nm);
  vector<Type> logrc(nm);
  vector<Type> rold(nm);
  vector<Type> logrold(nm);
  //vector<Type> rp(nm);
  //vector<Type> logrp(nm);
  for(int i=0; i<nm; i++){
    rold(i) =  CppAD::abs(gamma * m(i) / K);
    logrold(i) = log(rold(i));
    rc(i) = CppAD::abs(2.0 * rold(i) * (n - 1.0) / n);
    logrc(i) = log(rc(i));
    //rp(i) = abs(r(i) * (n - 1.0));
    //logrp(i) = log(rp(i));
    r(i) = m(i)/K * pow(n,(n/(n-1.0))); //abs(r(i) * (n - 1.0));
    logr(i) = log(r(i));
    //std::cout << " -- n: " << n << " -- gamma: " << gamma << n << " -- m(i): " << m(i)<< n << " -- K: " << K << " -- r(i): " << r(i) << " -- logr(i): " << logr(i) << std::endl;
  }
  for(int i=0; i<ns; i++){
    logrre(i) = log(mvec(i)/K * pow(n,(n/(n-1.0))));
  }
  Type BmsyB0 = pow(Type(1)/n,Type(1)/(n-Type(1)) );
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
    std::cout << "logobsC.size(): " << logobsC.size() << "  Cpred.size(): " << Cpred.size() << "  logobsI.size(): " << logobsI.size() << "  dt.size(): " << dt.size() << "  logF.size(): " << logF.size() << "  logu.rows(): " << logu.rows() << "  logu.cols(): " << logu.cols() << "  B.size(): " << B.size() << "  P.size(): " << P.size() << "  mvec.size(): " << mvec.size() << "  iq.size(): " << iq.size() << "  ic.size(): " << ic.size() << "  ir.size(): " << ir.size() << "  logsdf.size(): " << logsdf.size() << "  logFmsy.size(): " << logFmsy.size() << "  logFmsyvec.size(): " << logFmsyvec.size() << "  logBmsy.size(): " << logBmsy.size() << "  logBmsyvec.size(): " << logBmsyvec.size() << "  m.size(): " << m.size() << "  logphi.size(): " << logphi.size() << "  logphipar.size(): " << logphipar.size() << std::endl;
    std::cout << "Bmsys: " << Bmsys << std::endl;
    std::cout << "Fmsys: " << Fmsys << std::endl;
    std::cout << "MSYs: " << MSYs << std::endl;
    std::cout << "Bmsyd: " << Bmsyd << std::endl;
    std::cout << "Fmsyd: " << Fmsyd << std::endl;
    std::cout << "MSYd: " << MSYd << std::endl;
    std::cout << "isdf: " << isdf << std::endl;
  }
  // Calculate mvec if multiple rs are used (rarely the case).
  for(int i=0; i<ns; i++){
    ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
    if(dbg>1){
      std::cout << "-- i: " << i << " -- ind: " << ind << " -   mvec(i): " << mvec(i) << std::endl;
    }
    //mvec(i) = m(ind);
    //logFmsyvec(i) = logFmsy(ind);
    //logBmsyvec(i) = logBmsy(ind);
  }

  // Put wide smooth distributions on difficult parameters to stabilise optimisation.
  // Note that this contributes to the objective function, which therefore cannot be
  // regarded as a likelihood to be compared with likelihoods from other models.
  // Only apply these if there is no "manual" prior on the parameter and if stabilise == 1
  if (stabilise == 1){
    if (priorbkfrac(2) != 1){
      ans -= dnorm(logB(0) - logK, Type(-0.2234), Type(10.0), 1);
    }
    //ans -= dnorm(logB(0), Type(10.0), Type(10.0), 1);
    if(priorF(2) != 1){
      ans -= dnorm(logF(0), Type(-0.2234), Type(10.0), 1);
    }
    if(priorn(2) != 1){
      ans -= dnorm(logn, Type(0.6931472), Type(10.0), 1); // log(2) = 0.6931472
    }
    if(priorbeta(2) != 1){
      ans -= dnorm(logbeta, Type(0.0), Type(10.0), 1);
    }
    if(prioralpha(2) != 1){
      for(int i=0; i<nsdi; i++){
        ans -= dnorm(logalpha(i), Type(0.0), Type(10.0), 1);
      }
    }
    if(priorsde(2) != 1){
      ans -= dnorm(logsde, Type(-0.9162907), Type(10.0), 1); // log(0.4) = -0.9162907
    }
  }

  // PRIORS
  if(dbg > 0){
    std::cout << "PRIOR: priorn(0): " << priorn(0) << " -- priorn(1): " << priorn(1) << " -- priorn(2): " << priorn(2) << std::endl;
    std::cout << "PRIOR: priorr(0): " << priorr(0) << " -- priorr(1): " << priorr(1) << " -- priorr(2): " << priorr(2) << std::endl;
    std::cout << "PRIOR: priorK(0): " << priorK(0) << " -- priorK(1): " << priorK(1) << " -- priorK(2): " << priorK(2) << std::endl;
    std::cout << "PRIOR: priorm(0): " << priorm(0) << " -- priorm(1): " << priorm(1) << " -- priorm(2): " << priorm(2) << std::endl;
    std::cout << "PRIOR: priorq(0): " << priorq(0) << " -- priorq(1): " << priorq(1) << " -- priorq(2): " << priorq(2) << std::endl;
  }
  // Gamma priors
  // TMB uses shape and scale parameterisation, but spict uses shape and rate similar to jags.
  if(priorngamma(2)==1){
    ans-= dgamma(logn, priorngamma(0), 1.0/priorngamma(1), 1);
  }
  // Log-normal priors
  if(priorn(2) == 1){
    ans -= dnorm(logn, priorn(0), priorn(1), 1);
  }
  if(prioralpha(2) == 1){
    for(int i=0; i<nsdi; i++){
      ans -= dnorm(logalpha(i), prioralpha(0), prioralpha(1), 1);
    }
  }
  if(priorbeta(2) == 1){
    ans -= dnorm(logbeta, priorbeta(0), priorbeta(1), 1);
  }
  if((priorr(2) == 1) & (nm == 1)){
    ans-= dnorm(logr(0), priorr(0), priorr(1), 1); // Prior for logr
  }
  if(priorK(2) == 1){
    ans-= dnorm(logK, priorK(0), priorK(1), 1); // Prior for logK
  }
  if((priorm(2) == 1) & (nm == 1)){
    ans-= dnorm(logm(0), priorm(0), priorm(1), 1); // Prior for logm
  }
  if((priormu(2) == 1) & (nm == 1)){
    ans-= dnorm(mu, priormu(0), priormu(1), 1); // Prior for mu
  }
  for(int i=0; i<nq; i++){
    if(priorq(i,2) == 1){
      ans-= dnorm(logq(i), priorq(i,0), priorq(i, 1), 1); // Prior for logq - log-normal
    }
  }
  if(priorqf(2) == 1){
    ans-= dnorm(logqf, priorqf(0), priorqf(1), 1); // Prior for logqf
  }
  if(priorbkfrac(2) == 1){
    ans-= dnorm(logB(0) - logK, priorbkfrac(0), priorbkfrac(1), 1); // Prior for logbkfrac
  }
  if(priorsdb(2) == 1){
    ans-= dnorm(logsdb, priorsdb(0), priorsdb(1), 1); // Prior for logsdb
  }
  if(priorsdm(2) == 1){
    ans-= dnorm(logsdm, priorsdm(0), priorsdm(1), 1); // Prior for logsdm
  }
  if(priorsdf(2) == 1){
    for(int i=0; i<nsdf; i++){
      ans-= dnorm(logsdf(i), priorsdf(0), priorsdf(1), 1); // Prior for logsdf
    }
  }
  for(int i=0; i<nsdi; i++){
    if(priorsdi(i, 2) == 1){
      ans-= dnorm(logsdi(i), priorsdi(i, 0), priorsdi(i, 1), 1);  // Prior for logsdi
    }
  }
  if(priorsde(2) == 1){
    ans-= dnorm(logsde, priorsde(0), priorsde(1), 1); // Prior for logsde
  }
  if(priorsdc(2) == 1){
    ans-= dnorm(logsdc, priorsdc(0), priorsdc(1), 1); // Prior for logsdc
  }
  if(priorpsi(2) == 1){
    ans-= dnorm(logpsi, priorpsi(0), priorpsi(1), 1); // Prior for logsdm
  }
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
  if(priorFFmsy(2) == 1){
    ind = CppAD::Integer(priorFFmsy(4)-1);
    ans-= dnorm(logF(ind) - logFmsyvec(ind), priorFFmsy(0), priorFFmsy(1), 1); // Prior for logFFmsy
  }
  if(priorBmsyB0(2) == 1){
    ans -= dnorm( BmsyB0, priorBmsyB0(0),priorBmsyB0(1), 1); // prior for Bmsy/B0
  }

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

  int stepsPrYear = CppAD::Integer(1/dt[0]);
  vector<Type> SARphivec(nseasons);
  SARphivec.setZero();
  SARphivec(nseasons-1) = SARphi;

  using namespace density;
  ARk_t<Type> nldens(SARphivec);
  if(seasontype==3){
    ans += SCALE(nldens, sdSAR)(vector<Type>(SARvec));
    SIMULATE{
      if(simRandomEffects == 1) nldens.simulate(SARvec);
      REPORT(SARvec);
    }
  }
  // std::cout << "-- sdf2: " << sdf2 << std::endl;
  // std::cout << "-- sdf: " << sdf << std::endl;

  //vector<Type> logFs = logF
  vector<Type> logS(ns);
  vector<Type> logFpred(ns);
  if(simple==0){
    if(dbg>0){
      std::cout << "--- DEBUG: F loop start --- ans: " << ans << std::endl;
    }
    // Diffusion component of F
    int iisdf;
    for(int i=1; i<ns; i++){
      Type Fpredtmp = 0.0;
      iisdf = CppAD::Integer(isdf(i)) - 1;
      if (efforttype == 1){
        Fpredtmp = predictF1(logF(i-1), dt(i), sdf2(iisdf), delta, logeta);
      }
      if (efforttype == 2){
        Fpredtmp = predictF2(logF(i-1), dt(i), sdf2(iisdf), delta, logeta);
      }
      logFpred(i) = log( ffacvec(i) * Fpredtmp + fconvec(i) );
      likval = dnorm(logF(i), logFpred(i), sqrt(dt(i-1)) * sdf(iisdf), 1);
      ans-=likval;
      SIMULATE{
        if(simRandomEffects == 1) logF(i) = rnorm(logFpred(i), sqrt(dt(i-1)) * sdf(iisdf));
      }

      // DEBUGGING
      if(dbg>1){
        std::cout << "-- i: " << i << " -   logF(i-1): " << logF(i-1) << "  logF(i): " << logF(i) << "  ffacvec(i): " << ffacvec(i) << "  fconvec(i): " << fconvec(i) << "  iisdf: " << iisdf << "  sdf: " << sdf(iisdf) << "  efforttype: " << efforttype << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
    SIMULATE{
      REPORT(logF);
      vector<Type> trueF = exp(logFpred);
      REPORT(trueF);
    }

    // Seasonal component
    if(dbg>0){ std::cout << "-- seasontype: " << seasontype << std::endl; }
    for(int i=0; i<ns; i++) logS(i) = 0.0; // Initialise
    if(seasontype == 1 || seasontype == 3 ){
      // Spline
      int ind2, ind3;
      for(int i=0; i<ns; i++){
        ind2 = CppAD::Integer(seasonindex(i));
        ind3 = CppAD::Integer(seasonindex2(i));
        //logFs(i) += seasonspline(ind2);
        logS(i) += seasonspline(ind2);

        if(seasontype == 3) logS(i) += SARvec(ind3-1);
        // DEBUGGING
        if(dbg>1){
          //std::cout << "-- i: " << i << " -   logF(i): " << logF(i) << " logFs(i): " << logFs(i) << " ind2: " << ind2 << " seasonspline(ind2): " << seasonspline(ind2) << std::endl;
          std::cout << "-- i: " << i << " -   logF(i): " << logF(i) << " logS(i): " << logS(i) << " ind2: " << ind2 << " seasonspline(ind2): " << seasonspline(ind2) << std::endl;
        }
      }
    }
    if(seasontype == 2){
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
            SIMULATE{
              if(simRandomEffects == 1) logu(2*j+k, i) = rnorm(logupred(k), sduana);
            }
          }
          ans-=likval;
          // DEBUGGING
          if(dbg>0){
            std::cout << "-- i: " << i << " -   logu(0,i): " << logu(0,i) << "  logupred(0): " << logupred(0) << " -   logu(1,i): " << logu(1,i) << "  logupred(1): " << logupred(1) << "  sdu(j): " << sdu(j) << "  sduana: " << sduana << "  likval: " << likval << "  ans:" << ans << std::endl;
          }
        }
        //for(int i=0; i<ns; i++) logFs(i) += logu(2*j, i); // Sum diffusion and seasonal component
        for(int i=0; i<ns; i++) logS(i) += logu(2*j, i); // Sum diffusion and seasonal component
      }
    }
    SIMULATE{
      REPORT(logu);
    }

  } else {
    for(int i=0; i<ns; i++) logS(i) = -30; // If using simple set fishing mortality to something small.
  }
  vector<Type> F = exp(logS + logF); // This is the fishing mortality used to calculate catch
  vector<Type> logFs = log(F);

  SIMULATE{
    REPORT(logS);
    REPORT(logFs);
  }


  // GROWTH RATE (modelled as time-varying m)
  if (timevaryinggrowth == 1){
    if (dbg > 0){
      std::cout << "--- DEBUG: logmre loop start --- ans: " << ans << std::endl;
    }
    // Compare initial value with stationary distribution of OU
    likval = dnorm(logmre(0), Type(0.0), sdm/sqrt(2.0*psi), 1);
    SIMULATE{
      if(simRandomEffects == 1) logmre(0) = rnorm(Type(0.0), sdm / sqrt(2.0 * psi));
    }
    //likval = dnorm(logmre(0), logm(0), sdm/sqrt(2.0*psi), 1);
    ans -= likval;
    for (int i=1; i < ns; i++){
      Type logmrepred = predictm(logmre(i-1), dt(i-1), sdm2, psi);
      likval = dnorm(logmre(i), logmrepred, sqrt(dt(i-1))*sdm, 1);
      SIMULATE{
        if(simRandomEffects == 1) logmre(i) = rnorm(logmrepred, sqrt(dt(i-1)) * sdm);
      }
      //likval = dnorm(logmre(i), logmre(i-1), sqrt(dt(i-1))*sdm, 1);
      ans -= likval;
      // DEBUGGING
      if (dbg > 1){
        std::cout << "-- i: " << i << " -   logmre(i-1): " << logmre(i-1) << "  sdm: " << sdm << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
    SIMULATE{
      REPORT(logmre);
    }
  }

  // std::cout << "-- sdb2: " << sdb2 << std::endl;
  // std::cout << "-- sdb: " << sdb << std::endl;

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
      // Use naive approach
      Type Bpredtmp = exp(predictlogB(B(i), Ftmp, gamma, mvec(i), K, dt(i), n, sdb2) + 0.5*sdb2*dt(i)) - exp(logobsC(i));
      if(Bpredtmp < 0) Bpredtmp = 1e-8; // Ugly ugly ugly hack to avoid taking log of negative
      logBpred(i+1) = log(Bpredtmp);
      logFs(i) = logobsC(i) - logB(i); // Calculate fishing mortality
    }
    likval = dnorm(logBpred(i+1), logB(i+1), sqrt(dt(i))*sdb, 1);
    SIMULATE{
      if(simRandomEffects == 1){
        logB(i+1) = rnorm(logBpred(i+1), sqrt(dt(i)) * sdb);
        B(i+1) = exp(logB(i+1));
      }
    }
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logB(i+1): " << logB(i+1) << "  logBpred(i+1): " << logBpred(i+1) << "  sdb: " << sdb << "  likval: " << likval << "  ans:" << ans << std::endl;
    }
  }
  SIMULATE{
    REPORT(logB);
    vector<Type> trueB = exp(logBpred);
    REPORT(trueB);
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

  // EFFORT PREDICTIONS
  vector<Type> Fcumpred(nobsE);
  vector<Type> logFcumpred(nobsE);
  if(simple==0){
    for(int i=0; i<nobsE; i++){
      Fcumpred(i) = 0.0;
      // Sum effort contributions from each sub interval
      for(int j=0; j<ne(i); j++){
        ind = CppAD::Integer(ie(i)-1) + j; // minus 1 because R starts at 1 and c++ at 0
        Fcumpred(i) += F(ind) * dt(ind);
      }
      logFcumpred(i) = log(Fcumpred(i));
      // DEBUGGING
      if(dbg>1){
        std::cout << "-- i: " << i << " -  ind: " << ind << "  logFcumpred(i): " << logFcumpred(i) << std::endl;
      }
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
      if(robflagc==1){
        likval = log(pp*dnorm(logCpred(i), logobsC(i), stdevfacc(i)*sdc, 0) + (1.0-pp)*dnorm(logCpred(i), logobsC(i), robfac*stdevfacc(i)*sdc, 0));
        SIMULATE{
          Type uu = runif(0.0,1.0);
          if(uu < pp){
            logobsC(i) = log(rnorm(logCpred(i), stdevfacc(i) * sdc));
          }else{
            logobsC(i) = log(rnorm(logCpred(i), robfac * stdevfacc(i) * sdc));
          }
        }

      } else {
        likval = dnorm(logCpred(i), logobsC(i), stdevfacc(i)*sdc, 1);
        SIMULATE{
          logobsC(i) = rnorm(logCpred(i), stdevfacc(i) * sdc);
        }
      }
      ans-= keep(inds) * likval;
      // DEBUGGING
      if(dbg>1){
        std::cout << "-- i: " << i << " -   logobsC(i): " << logobsC(i) << " -   stdevfacc(i): " << stdevfacc(i) << "  sdc: " << sdc << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
    SIMULATE{
      vector<Type> obsC = exp(logobsC);
      REPORT(obsC);
      vector<Type> trueC = exp(logCpred);
      REPORT(trueC);
    }
  }

  // EFFORT
  if(simple==0){
    if(dbg>0){
      std::cout << "--- DEBUG: Epred loop start --- ans: " << ans << std::endl;
    }
    for(int i=0; i<nobsE; i++){
      logEpred(i) = logFcumpred(i) - logqf; // E = 1/q * integral{F_t dt}
      if(robflage==1){
        likval = log(pp*dnorm(logEpred(i), logobsE(i), stdevface(i)*sde, 0) + (1.0-pp)*dnorm(logEpred(i), logobsE(i), robfac*stdevface(i)*sde, 0));
        SIMULATE{
          Type uu = runif(0.0, 1.0);
          if(uu < pp){
            logobsE(i) = log(rnorm(logEpred(i), stdevface(i) * sde));
          }else{
            logobsE(i) = log(rnorm(logEpred(i), robfac * stdevface(i) * sde));  // same robfac and pp for C and E?
          }
        }
      } else {
        likval = dnorm(logEpred(i), logobsE(i), stdevface(i)*sde, 1);
        SIMULATE{
          logobsE(i) = rnorm(logEpred(i), stdevface(i) * sde);
        }
      }
      inds = CppAD::Integer(ise(i)-1);
      ans-= keep(inds) * likval;
      // DEBUGGING
      if(dbg>1){
        std::cout << "-- i: " << i << " -   logobsE(i): " << logobsE(i) << " -   stdevface(i): " << stdevface(i) << "  sde: " << sde << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
    SIMULATE{
      vector<Type> obsE = exp(logobsE);
      REPORT(obsE);
      vector<Type> trueE = exp(logEpred);
      REPORT(trueE);
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
    if(robflagi(indsdi)==1){
      likval = log(pp*dnorm(logobsI(i), logIpred(i), stdevfaci(i)*sdi(indsdi), 0) + (1.0-pp)*dnorm(logobsI(i), logIpred(i), robfac*stdevfaci(i)*sdi(indsdi), 0));
      SIMULATE{
        Type uu = runif(0.0,1.0);
        if(uu < pp){
          logobsI(i) = rnorm(logIpred(i), stdevfaci(i) * sdi(indsdi));
        }else{
          logobsI(i) = rnorm(logIpred(i), robfac * stdevfaci(i) * sdi(indsdi));
        }
      }
    } else {
      likval = dnorm(logobsI(i), logIpred(i), stdevfaci(i)*sdi(indsdi), 1);
      SIMULATE{
        logobsI(i) = rnorm(logIpred(i), stdevfaci(i) * sdi(indsdi));
      }
    }
    ans-= keep(inds) * likval * iuse(i);
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -  ind: " << ind << " -  indq: " << indq << " -  indsdi: " << indsdi << " -  inds: " << inds << " -   logobsI(i): " << logobsI(i) << "  logIpred(i): " << logIpred(i) << "  stdevfaci(i): " << stdevfaci(i) << "  likval: " << likval << "  sdi: " << sdi << "  ans:" << ans << std::endl;
    }
  }
  SIMULATE{
    vector<Type> obsI = exp(logobsI);
    REPORT(obsI);
    vector<Type> trueI = exp(logIpred);
    REPORT(trueI);
  }

  SIMULATE{
    // Put obs back into obssrt
    for(int i=0; i<nobsC; i++){
      ind = CppAD::Integer(isc(i)-1);
      obssrt(ind) = logobsC(i);
    }
    for(int i=0; i<nobsI; i++){
      ind = CppAD::Integer(isi(i)-1);
      obssrt(ind) = logobsI(i);
    }
    for(int i=0; i<nobsE; i++){
      ind = CppAD::Integer(ise(i)-1);
      obssrt(ind) = logobsE(i);
    }
    REPORT(obssrt);
  }


  /*
    --- ONE-STEP-AHEAD PREDICTIONS ---
  */

  if(dbg > 0){
    std::cout << "--- DEBUG: ONE-STEP-AHEAD CATCH PREDICTIONS --- ans: " << ans << std::endl;
    std::cout << "-- dtpredcnsteps: " << dtpredcnsteps << "  dtpredcinds.size(): " << dtpredcinds.size() <<std::endl;

  }
  // Catch prediction
  Type Cp = 0.0;
  for(int i=0; i<dtpredcnsteps; i++){
    ind = CppAD::Integer(dtpredcinds(i)-1);
    if(dbg>1){
      std::cout << "-- i: " << i << " -  dtpredcinds(i)-1: " << ind << std::endl;
    }
    Cp += Cpredsub(ind);
  }
  Type logCp = log(Cp);

  if(dbg > 0){
    std::cout << "--- DEBUG: ONE-STEP-AHEAD EFFORT PREDICTIONS --- ans: " << ans << std::endl;
    std::cout << "-- dtpredensteps: " << dtpredensteps << "  dtpredeinds.size(): " << dtpredeinds.size() <<std::endl;
  }

  // Effort prediction
  Type Fp = 0.0;
  for(int i=0; i<dtpredensteps; i++){
    ind = CppAD::Integer(dtpredeinds(i)-1);
    if(dbg>1){
      std::cout << "-- i: " << i << " -  dtpredeinds(i)-1: " << ind << std::endl;
    }
    Fp += F(ind);
  }
  Type logEp = log(Fp) - logqf;
  Type Ep = exp(logEp);

  if(dbg > 0){
    std::cout << "--- DEBUG: Biomass and F at the end of the prediction time interval --- ans: " << ans << std::endl;
  }

  // Biomass and F at the end of the prediction time interval
  int pind = CppAD::Integer(dtprediind-1);
  if(dbg > 0){
    std::cout << "-- pind: " << pind << std::endl;
  }
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

  if(dbg > 0){
    std::cout << "--- DEBUG: Biomass and fishing mortality at last time point --- ans: " << ans << std::endl;
  }

  // Biomass and F at the start of the prediction time interval
  int mind = CppAD::Integer(dtpredcinds(0)-1);
  if(dbg > 3){
    std::cout << "-- mind: " << mind << std::endl;
  }
  Type Bm = B(mind);
  Type logBm = log(Bm);
  Type logBmBmsy = logBm - logBmsyvec(mind);
  Type logBmK = logBm - logK;
  Type logFm = logFs(mind);
  Type logFmFmsy = logFm - logFmsyvec(mind);

  // Biomass and fishing mortality at last time point
  Type logBl = logB(indlastobs-1);
  Type logBlBmsy = logBl - logBmsyvec(indlastobs-1);
  Type logBlK = logBl - logK;
  Type logFl = logFs(indlastobs-1);
  Type logFlFmsy = logFl - logFmsyvec(indlastobs-1);

  if(dbg > 0){
    std::cout << "--- DEBUG: Calculate relative levels of biomass and fishing mortality --- ans: " << ans << std::endl;
  }

  // Calculate relative levels of biomass and fishing mortality
  vector<Type> logBBmsy(ns);
  vector<Type> logFFmsy(ns);
  Type meanB = exp(logB).sum()/ns;
  Type meanF = exp(logF).sum()/ns;
  vector<Type> logBrel(ns);
  vector<Type> logFrel(ns);

  for(int i=0; i<ns; i++){
    logBBmsy(i) = logB(i) - logBmsyvec(i);
    logFFmsy(i) = logFs(i) - logFmsyvec(i);
    logBrel(i) = logB(i) - log(meanB);
    logFrel(i) = logF(i) - log(meanF);
  }

  //std::cout << "logFFmsy: " << logFFmsy << std::endl;

  //
  Type logbkfrac = logB(0) - logK;

  if(dbg > 0){
    std::cout << "--- DEBUG: Calculations done, doing ADREPORTS --- ans: " << ans << std::endl;
  }

  vector<Type> logFnotS(ns);
  vector<Type> logFFmsynotS(ns);

  Type meanS = (exp(logFs - logF)).sum() / logF.size();
  for(int i=0; i<ns; i++){
    logFnotS(i) = logF(i) + log(meanS);
    logFFmsynotS(i) = logFnotS(i) - logFmsyvec(i);
  }
  Type logFlFmsynotS = logFnotS(indlastobs-1) - logFmsyvec(indlastobs-1);
  Type logFpFmsynotS = logFnotS(pind) - logFmsyvec(pind);
  Type logFmnotS = logFnotS(mind);
  Type logFmFmsynotS = logFmnotS - logFmsyvec(mind);

  // Report the sum of reference points -- can be used to calculate their covariance without using ADreport with covariance.
  Type logBmsyPluslogFmsy = logBmsy(logBmsy.size()-1) + logFmsy(logFmsy.size()-1);

  // ADREPORTS
  if(reportmode == 0){
    ADREPORT(Bmsy);
    ADREPORT(Bmsyd);
    ADREPORT(Bmsys);
    ADREPORT(Bmsy2);
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
    ADREPORT(logEp);
    // PARAMETERS
    ADREPORT(r);
    ADREPORT(logr);
    ADREPORT(rc);
    ADREPORT(logrc);
    //ADREPORT(rp);
    //ADREPORT(logrp);
    ADREPORT(rold);
    ADREPORT(logrold);
    ADREPORT(K);
    ADREPORT(q);
    //ADREPORT(logq);
    ADREPORT(logq2);
    ADREPORT(p);
    ADREPORT(gamma);
    ADREPORT(m);
    ADREPORT(sdf);
    ADREPORT(sdc);
    ADREPORT(sde);
    ADREPORT(sdb);
    ADREPORT(sdi);
    ADREPORT(isdf2);
    ADREPORT(isdc2);
    ADREPORT(isde2);
    ADREPORT(isdb2);
    ADREPORT(isdi2);
    ADREPORT(logalpha);
    ADREPORT(logbeta);
    ADREPORT(BmsyB0);
    if(reportall){
      // These reports are derived from the random effects and are therefore vectors. TMB calculates the covariance of all sdreports leading to a very large covariance matrix which may cause memory problems.
      // B
      ADREPORT(logBBmsy);
      if(reportRel) ADREPORT(logBrel);
      // F
      ADREPORT(logFFmsy); // Vector of size ns
      ADREPORT(logFs);    // Vector of size ns
      if(reportRel) ADREPORT(logFrel);
      // C
      ADREPORT(logCpred);
      // I
      ADREPORT(logIpred);
      // E
      ADREPORT(logEpred);
      // Time varying growth
      if ((timevaryinggrowth == 1) | (logmcovflag == 1)){
        ADREPORT(logrre); // r random effect
        ADREPORT(logFmsyvec);
        ADREPORT(logMSYvec);
      }
      ADREPORT(logFnotS);
      ADREPORT(logFFmsynotS);
    }
    ADREPORT( logBmsyPluslogFmsy ) ;

    ADREPORT(Bm);
    ADREPORT(logBm);
    ADREPORT(logBmBmsy);
    ADREPORT(logBmK);
    ADREPORT(logFm);
    ADREPORT(logFmFmsy);
    ADREPORT(logFlFmsynotS);
    ADREPORT(logFmFmsynotS);
    ADREPORT(logFpFmsynotS);

  }else if(reportmode == 1){
    ADREPORT(logFm);
    ADREPORT(logFmsy);
    ADREPORT(logBmsy);
    ADREPORT(logFpFmsynotS);
    ADREPORT(logBpBmsy);
    ADREPORT(logFmFmsynotS);
    ADREPORT(logBmBmsy);
    ADREPORT(logCp);
  }else if(reportmode == 2){
    ADREPORT(logCp);
  }

  // REPORTS (these don't require sdreport to be output)
  REPORT(Cp);
  REPORT(Ep);
  REPORT(logIp);
  REPORT(MSY);
  REPORT(Bmsy);
  REPORT(Fmsy);
  // REPORT(stochmsy);   // otherwise checkConsistency doesn't work
  // REPORT(logBBmsy);
  // REPORT(logFFmsy);
  // REPORT(logB);
  // REPORT(logF);

  return ans;
}
