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
#include <TMB.hpp>

/* Predict biomass */
template<class Type>
Type predictlogB(const Type &B0, const Type &F, const Type &gamma, const Type &m, const Type &K, const Type &dt, const Type &n, const Type &sdb2)
{
  // Euler discretised Lamperti transformed Pella-Tomlinson surplus production model in Fletcher (1978) form.
  return log(B0) + (gamma*m/K - gamma*m/K*pow(B0/K, n-1.0) - F - 0.5*sdb2)*dt;
}

/* Predict F1 */
template<class Type>
Type predictE1(const Type &logE0, const Type &dt, const Type &sdf2)
{
  return exp(logE0);
}

/* Predict F2 */
template<class Type>
Type predictE2(const Type &logE0, const Type &dt, const Type &sdf2)
{
  return exp(logE0 - 0.5*sdf2*dt);
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
  DATA_VECTOR(index2sdb);    // 
  DATA_VECTOR(index2sdi);    // 
  DATA_VECTOR(index2q);    // 
  DATA_INTEGER(nfleets);       // 
  DATA_VECTOR(ic);             // Vector such that B(ic(i)) is the state at the start of obsC(i)
  DATA_VECTOR(nc);             // nc(i) gives the number of time intervals obsC(i) spans
  DATA_VECTOR(icpred);         // 
  DATA_VECTOR(iff);            // 
  DATA_VECTOR(ie);             // Vector such that E(ie(i)) is the state at the start of obsE(i)
  DATA_VECTOR(ifleet);         // length nobsE
  DATA_VECTOR(ne);             // ne(i) gives the number of time intervals obsE(i) spans
  DATA_VECTOR(ii);             // A vector such that B(ii(i)) is the state corresponding to I(i)
  DATA_VECTOR(ib);             // 
  DATA_VECTOR(iq);             // A vector such that iq(i) is the index number corresponding to I_iq(i)
  DATA_VECTOR(isdi);           // A vector such that isdi(i) is the index number corresponding to I_isdi(i)
  DATA_VECTOR(isdc);           // 
  DATA_VECTOR(isde);           // 
  //DATA_VECTOR(ir);             // A vector indicating when the different rs should be used
  DATA_VECTOR(seasons);        // A vector of length ns indicating to which season a state belongs
  DATA_VECTOR(seasonindex);    // A vector of length ns giving the number stepped within the current year
  DATA_MATRIX(splinemat);      // Design matrix for the seasonal spline
  DATA_MATRIX(splinematfine);  // Design matrix for the seasonal spline on a fine time scale to get spline uncertainty
  DATA_SCALAR(omega);          // Period time of seasonal SDEs (2*pi = 1 year period)
  DATA_SCALAR(seasontype);     // Variable indicating whether to use 1=spline, 2=coupled SDEs
  DATA_SCALAR(efforttype);     // Variable indicating whether to use 1=spline, 2=coupled SDEs
  DATA_VECTOR(ffacvec);        // Management factor each year multiply the predicted F with ffac
  DATA_VECTOR(fconvec);        // Management factor each year add this constant to the predicted F
  DATA_VECTOR(indpred);        // A vector indicating when the management factor should be applied
  DATA_MATRIX(targetmap);      // Matrix where nrows = nqf, first column is stock index, second column is fleet index
  DATA_SCALAR(robflagc);       // If 1 use robust observation error for catches
  DATA_SCALAR(robflagi);       // If 1 use robust observation error for index
  DATA_SCALAR(robflage);       // If 1 use robust observation error for effort
  DATA_INTEGER(stochmsy);      // Use stochastic msy?
  //DATA_SCALAR(effortflag);     // If effortflag == 1 use effort data, else use index data

  // Priors
  DATA_VECTOR(priorn);         // Prior vec for n, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorr);         // Prior vec for r, [log(mean), stdev in log, useflag]
  //DATA_VECTOR(priorrp);        // Prior vec for rp, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorK);         // Prior vec for K, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorm);         // Prior vec for m, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorq);         // Prior vec for q, [log(mean), stdev in log, useflag]
  DATA_VECTOR(prioriqgamma);   // Prior vec for q, inverse gamma distribution, [shape, rate, useflag]
  DATA_VECTOR(priorqf);        // Prior vec for qf, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorbkfrac);    // Prior vec for B0/K, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorsdb);       // Prior vec for sdb, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorisdb2gamma);// Prior vec for sdb2, inv. gamma distribution, [shape, rate, useflag]
  DATA_VECTOR(priorsdf);       // Prior vec for sdf, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorisdf2gamma);// Prior vec for sdf2, inv. gamma distribution, [shape, rate, useflag]
  DATA_VECTOR(priorsdi);       // Prior vec for sdi, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorisdi2gamma);// Prior vec for sdi2, inv. gamma distribution, [shape, rate, useflag]
  DATA_VECTOR(priorsde);       // Prior vec for sde, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorisde2gamma);// Prior vec for sde2, inv. gamma distribution, [shape, rate, useflag]
  DATA_VECTOR(priorsdc);       // Prior vec for sdc, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorisdc2gamma);// Prior vec for sdc2, inv. gamma distribution, [shape, rate, useflag]
  DATA_VECTOR(prioralpha);     // Prior vec for alpha, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorbeta);      // Prior vec for beta, [log(mean), stdev in log, useflag]
  DATA_VECTOR(priorB);         // Prior vec for B, [log(mean), stdev in log, useflag, year, ib]
  DATA_VECTOR(priorF);         // Prior vec for F, [log(mean), stdev in log, useflag, year, if]
  DATA_VECTOR(priorBBmsy);     // Prior vec for B/Bmsy, [log(mean), stdev in log, useflag, year, ib]
  DATA_VECTOR(priorFFmsy);     // Prior vec for F/Fmsy, [log(mean), stdev in log, useflag, year, if]

  // Options
  DATA_SCALAR(simple);         // If simple=1 use simple model (catch assumed known, no F process)
  DATA_SCALAR(dbg);            // Debug flag, if == 1 then print stuff.

  // PARAMETERS
  PARAMETER_VECTOR(logm);      // MSY, length nstocks
  PARAMETER_VECTOR(logK);      // Carrying capacity, length nstocks
  PARAMETER_VECTOR(logq);      // Catchability for index, length sum(nindex)
  PARAMETER_VECTOR(logqf);     // Catchability for effort, length nqf
  PARAMETER_VECTOR(logn);      // Pella-Tomlinson exponent, length nstocks
  PARAMETER_VECTOR(logsdb);    // Stdev in B process, length nstocks
  PARAMETER_VECTOR(logsdu);    // Stdev in seasonal component of F process, length nqf
  PARAMETER_VECTOR(logsdf);    // Stdev in diffusion component of F process, length nqf
  PARAMETER_VECTOR(logsdi);    // Stdev of index obs noise, length sum(nindex)
  PARAMETER_VECTOR(logsde);    // Stdev of effort obs noise, length nqf
  PARAMETER_VECTOR(logsdc);    // Stdev of catch obs noise, length nqf
  PARAMETER_VECTOR(logphi);    // Season levels of F when using spline based seasons
  PARAMETER_VECTOR(loglambda); // Damping variable when using seasonal SDEs
  PARAMETER(logitpp);          // Proportion of narrow distribution when using robust obs err.
  PARAMETER(logp1robfac);      // Coefficient to the standard deviation of robust observation error distribution
  PARAMETER_MATRIX(logE);      // Diffusion component of F in log
  PARAMETER_MATRIX(logu);      // Seasonal component of F in log
  PARAMETER_MATRIX(logB);      // Biomass in log

  //std::cout << "expmosc: " << expmosc(lambda, omega, 0.1) << std::endl;
   if (dbg > 0){
     std::cout << "==== DATA read, now calculating derived quantities ====" << std::endl;
   }

  int ind;
  int find;
  int sind;
  int bind;
  int iind;
  // Distribute sorted observations into logobsC and logobsI vectors
  //int nobsC = isc.size();
  vector<Type> logobsC(nobsC);
  for (int i=0; i < nobsC; i++){ 
    ind = CppAD::Integer(isc(i)-1);
    logobsC(i) = obssrt(ind); 
  }
  //int nobsI = isi.size();
  vector<Type> logobsI(nobsI);
  for (int i=0; i < nobsI; i++){ 
    ind = CppAD::Integer(isi(i)-1);
    logobsI(i) = obssrt(ind);
  }
  vector<Type> logobsE(nobsE);
  for (int i=0; i < nobsE; i++){ 
    ind = CppAD::Integer(ise(i)-1);
    logobsE(i) = obssrt(ind);
  }

  // Length of vectors
  int nstocks = logm.size();
  int nq = logq.size();
  int nsdu = logsdu.size();
  int nsdi = logsdi.size();
  int nsde = logsde.size();
  int nobsCp = ic.size();
  int ns = logE.cols(); // is this allowed?
  int nqf = logqf.size();
  int nindex = index2sdb.size();
  //int nfleets = ifleet.max(); // is this allowed?
  //int nfleets = max(targetmap.col(1)); // is this allowed?

  // Create fishing mortality, logF
  matrix<Type> logF(nqf, ns);
  matrix<Type> F(nqf, ns);
  for (int k=0; k < nqf; k++){
    find = CppAD::Integer(targetmap(k, 1) - 1); // minus 1 because R starts at 1 and c++ at 0
    for (int i=0; i < ns; i++){
      logF(k, i) = logqf(k) + logE(find, i);
    }
  }
  for (int k=0; k < nqf; k++){
    for (int i=0; i < ns; i++){
      F(k, i) = exp(logF(k, i));
    }
  }

  // Transform parameters
  vector<Type> logphipar(logphi.size()+1);
  logphipar(0) = 0.0; // The first logphi is set to 0, the rest are estimated relative to this.
  for (int i=1; i < logphipar.size(); i++){ 
    logphipar(i) = logphi(i-1); 
  }
  vector<Type> temp = splinemat.col(0);
  vector<Type> seasonspline(temp.size());
  seasonspline = splinemat * logphipar;
  vector<Type> tempfine = splinematfine.col(0);
  vector<Type> seasonsplinefine(tempfine.size());
  seasonsplinefine = splinematfine * logphipar;
  Type pp = 1.0 / (1.0 + exp(-logitpp));
  Type robfac = 1.0 + exp(logp1robfac);
  vector<Type> m(nstocks); 
  for (int i=0; i < nstocks; i++){ 
    m(i) = exp(logm(i)); 
  }
  vector<Type> K(nstocks); 
  for (int i=0; i < nstocks; i++){ 
    K(i) = exp(logK(i)); 
  }
  //Type K = exp(logK);
  vector<Type> q(nq);
  for (int i=0; i < nq; i++){ 
    q(i) = exp(logq(i)); 
  }

  vector<Type> logq2(nq); // This is to be able to compare with Polacheck et al. (1993)
  for (int i=0; i < nq; i++){ 
    logq2(i) = log(100) + logq(i); 
  }
  vector<Type> qf(nqf);
  for (int i=0; i < nqf; i++){ 
    qf(i) = exp(logqf(i)); 
  }
  //Type qf = exp(logqf);
  vector<Type> n(nstocks);
  for (int i=0; i < nstocks; i++){ 
    n(i) = exp(logn(i)); 
  }
  //Type n = exp(logn);
  vector<Type> gamma(nstocks);
  for (int i=0; i < nstocks; i++){ 
    gamma(i) = pow(n(i), n(i)/(n(i)-1.0)) / (n(i)-1.0);
  }
  //Type gamma = pow(n, n/(n-1.0)) / (n-1.0);
  vector<Type> lambda(nqf);
  for (int i=0; i < nqf; i++){ 
    lambda(i) = exp(loglambda(i)); 
  }
  //Type lambda = exp(loglambda);
  vector<Type> sdf(nqf);
  vector<Type> sdf2(nqf);
  vector<Type> isdf2(nqf);
  for (int i=0; i < nqf; i++){ 
    sdf(i) = exp(logsdf(i)); 
    sdf2(i) = sdf(i) * sdf(i);
    isdf2(i) = 1.0 / sdf2(i);
  }
  // sdu
  vector<Type> sdu(nqf);
  for (int i=0; i < nqf; i++){ 
    sdu(i) = exp(logsdu(i)); 
  }
  // sdb
  vector<Type> sdb(nstocks);
  vector<Type> sdb2(nstocks);
  vector<Type> isdb2(nstocks);
  for (int i=0; i < nstocks; i++){ 
    sdb(i) = exp(logsdb(i)); 
    sdb2(i) = sdb(i) * sdb(i);
    isdb2(i) = 1.0 / sdb2(i);
  }
  // sdi
  vector<Type> sdi(nsdi);
  vector<Type> sdi2(nsdi);
  vector<Type> isdi2(nsdi);
  for (int i=0; i < nsdi; i++){
    sdi(i) = exp(logsdi(i)); 
    sdi2(i) = sdi(i) * sdi(i);
    isdi2(i) = 1.0 / sdi2(i);
  }
  // alpha
  vector<Type> alpha(nindex);
  for (int i=0; i < nindex; i++){
    bind = CppAD::Integer(index2sdb(i)-1);
    iind = CppAD::Integer(index2sdi(i)-1);
    alpha(i) = sdi(iind) / sdb(bind);
  }
  //vector<Type> alpha = sdi/sdb;
  vector<Type> logalpha = log(alpha);
  // sde
  vector<Type> sde(nsde);
  vector<Type> sde2(nsde);
  vector<Type> isde2(nsde);
  for (int i=0; i < nsde; i++){
    sde(i) = exp(logsde(i)); 
    sde2(i) = sde(i) * sde(i);
    isde2(i) = 1.0 / sde2(i);
  }
  // sdc
  vector<Type> sdc(nqf);
  vector<Type> sdc2(nqf);
  vector<Type> isdc2(nqf);
  for (int i=0; i < nqf; i++){
    sdc(i) = exp(logsdc(i)); 
    sdc2(i) = sdc(i) * sdc(i);
    isdc2(i) = 1.0 / sdc2(i);
  }
  // beta
  vector<Type> beta(nqf);
  for (int i=0; i < nqf; i++){
    beta(i) = sdc(i) / sdf(i);
  }
  vector<Type> logbeta = log(beta);

  // Put wide smooth distributions on difficult parameters to stabilise optimisation.
  // Note that this contributes to the objective function, which therefore cannot be 
  // regarded as a likelihood to be compared with likelihoods from other models.
  for (int i=0; i < nstocks; i++){
    ans -= dnorm(logB(i, 0) - logK(i), Type(-0.2234), Type(10.0), 1);
    ans -= dnorm(logn(i), Type(0.6931472), Type(10.0), 1); // log(2) = 0.6931472
  }
  for (int i=0; i < nqf; i++){
    ans -= dnorm(logF(i, 0), Type(-0.2234), Type(10.0), 1);
    ans -= dnorm(logbeta(i), Type(0.0), Type(10.0), 1);
  }
  for (int i=0; i < nindex; i++){
    ans -= dnorm(logalpha(i), Type(0.0), Type(10.0), 1); 
  }
  for (int i=0; i < nsde; i++){
    ans -= dnorm(logsde(i), Type(-0.9162907), Type(10.0), 1); // log(0.4) = -0.9162907
  }

  // Initialise vectors
  matrix<Type> B(nstocks, ns); // = exp(logB);
  for (int si=0; si < nstocks; si++){
    for (int i=0; i < ns; i++){
      B(si, i) = exp(logB(si, i));
    }
  }
  vector<Type> mvec(ns);
  matrix<Type> logBmsyvec(nstocks, ns);
  matrix<Type> logFmsyvec(nstocks, ns);
  vector<Type> Cpred(nobsCp);
  for (int i=0; i < nobsCp; i++){
    Cpred(i) = 0.0; 
  }
  vector<Type> logCpred(nobsCp);
  vector<Type> logIpred(nobsI);
  vector<Type> logEpred(nobsE);

  // Reference points
  vector<Type> p(nstocks);
  for (int si=0; si < nstocks; si++){
    p(si) = n(si) - 1.0; // Is this allowed?
  }
  vector<Type> Bmsyd(nstocks);
  vector<Type> Fmsyd(nstocks);
  vector<Type> MSYd = m;
  vector<Type> Bmsys(nstocks);
  vector<Type> Fmsys(nstocks);
  vector<Type> MSYs(nstocks);
  int flag;
  for (int si=0; si < nstocks; si++){
    flag = asDouble(n(si)) > 1; // Cast n as double to calc flag
    // Deterministic reference points
    Bmsyd(si) = K(si) * pow(1.0/n(si), 1.0/(n(si)-1.0));
    Fmsyd(si) = MSYd(si) / Bmsyd(si);
    // Stochastic reference points (NOTE: only proved for n > 1, Bordet and Rivest (2014))
    // The stepfun ensures that stochastic reference points are only used if n > 1.
    Type BmsyStochContr = (1.0 - (1.0 + Fmsyd(si)*(p(si)-1.0)/2.0)*sdb2(si) / (Fmsyd(si)*pow(2.0-Fmsyd(si), 2.0)));
    //Bmsys(si) = Bmsyd(si) * pow(BmsyStochContr, stepfun(p));
    //Bmsys(si) = Bmsyd(si) * (1.0 - (1.0 + Fmsyd(si)*(p-1.0)/2.0)*sdb2 / (Fmsyd(si)*pow(2.0-Fmsyd(si), 2.0)));
    Type FmsyStochContr = (p(si)*(1.0-Fmsyd(si))*sdb2(si)) / pow(2.0-Fmsyd(si), 2.0);
    //Fmsys(si) = Fmsyd(si) - stepfun(p) * FmsyStochContr;
    //Fmsys(si) = Fmsyd(si) - (p*(1.0-Fmsyd(si))*sdb2) / pow(2.0-Fmsyd(si), 2.0);
    Type MSYstochContr = (1.0 - ((p(si)+1.0)/2.0*sdb2(si)) / (1.0 - pow(1.0-Fmsyd(si), 2.0)));
    //MSYs(si) = MSYd(si) * pow(MSYstochContr, stepfun(p));
    //MSYs(si) = MSYd(si) * (1.0 - ((p+1.0)/2.0*sdb2) / (1.0 - pow(1.0-Fmsyd(si), 2.0)));


    //flag = asDouble(n) > 1 & asDouble(BmsyStochContr) > 0;
    if (flag){
      Bmsys(si) = Bmsyd(si) * BmsyStochContr;
      Fmsys(si) = Fmsyd(si) - FmsyStochContr;
      MSYs(si) = MSYd(si) * MSYstochContr;
    } else {
      Bmsys(si) = Bmsyd(si);
      Fmsys(si) = Fmsyd(si);
      MSYs(si) = MSYd(si);
    }
    //std::cout << "flag: " << flag << std::endl;
    //std::cout << "BmsyStochContr: " << BmsyStochContr << std::endl;
    //std::cout << "Bmsyd(si): " << Bmsyd(si) << std::endl;
    //std::cout << "Bmsys(si): " << Bmsys(si) << std::endl;
  }

  // log reference points
  vector<Type> logBmsyd = log(Bmsyd);
  vector<Type> logMSYd = log(MSYd);
  vector<Type> logFmsyd = log(Fmsyd);
  vector<Type> logBmsys = log(Bmsys);
  vector<Type> logFmsys = log(Fmsys);
  vector<Type> logMSYs = log(MSYs);
  // Used reference points
  vector<Type> Bmsy(nstocks);
  vector<Type> MSY(nstocks);
  vector<Type> Fmsy(nstocks);
  vector<Type> logBmsy(nstocks);
  vector<Type> logFmsy(nstocks);
  vector<Type> logMSY(nstocks);

  vector<Type> Bmsy2(nstocks);
  for (int si=0; si < nstocks; si++){
    flag = asDouble(n(si)) > 1; // Cast n as double to calc flag
    if (flag){
      Bmsy2(si) = Bmsys(si);
    } else {
      Bmsy2(si) = Bmsyd(si);
    }
  }

  if (stochmsy == 1){
    // Use stochastic reference points
    Bmsy = Bmsys;
    MSY = MSYs;
    Fmsy = Fmsys;
    logBmsy = logBmsys;
    logFmsy = logFmsys;
    logMSY = logMSYs;
    for (int si=0; si < nstocks; si++){
      for (int i=0; i < ns; i++){
	//ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
	logBmsyvec(si, i) = logBmsys(si);
	logFmsyvec(si, i) = logFmsys(si);
      }
    }
  } else {
    // Use deterministic reference points
    Bmsy = Bmsyd;
    MSY = MSYd;
    Fmsy = Fmsyd;
    logBmsy = logBmsyd;
    logFmsy = logFmsyd;
    logMSY = logMSYd;
    for (int si=0; si < nstocks; si++){
      for (int i=0; i < ns; i++){
	//ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
	logBmsyvec(si, i) = logBmsyd(si);
	logFmsyvec(si, i) = logFmsyd(si);
      }
    }
  }

  // These quantities are calculated to enable comparison with the Polacheck et al (1993) parameter estimates
  /*
  vector<Type> Emsy(nq);
  vector<Type> Emsy2(nq);
  for (int i=0; i<nq; i++){ 
    Emsy(i) = Fmsy(0) / q(i); 
    Emsy2(i) = Fmsy(0) / exp(logq2(i))*1.0e4; // Used for the results of the albacore data set
  }
  vector<Type> logEmsy = log(Emsy);
  vector<Type> logEmsy2 = log(Emsy2);
  */

  // Calculate growth rate
  //Type sign = 1.0;
  //if (n < 1.0) sign = -1.0; // Following Fletcher (1978)
  vector<Type> r(nstocks);
  vector<Type> logr(nstocks);
  vector<Type> rc(nstocks);
  vector<Type> logrc(nstocks);
  vector<Type> rold(nstocks);
  vector<Type> logrold(nstocks);
  //vector<Type> rp(nstocks);
  //vector<Type> logrp(nstocks);
  for (int i=0; i < nstocks; i++){ 
    rold(i) =  abs(gamma(i) * m(i) / K(i));
    logrold(i) = log(rold(i)); 
    rc(i) = abs(2.0 * rold(i) * (n(i) - 1.0) / n(i));
    logrc(i) = log(rc(i)); 
    //rp(i) = abs(r(i) * (n - 1.0));
    //logrp(i) = log(rp(i)); 
    r(i) = m(i) / K(i) * pow(n(i), (n(i)/(n(i)-1.0))); //abs(r(i) * (n - 1.0));
    logr(i) = log(r(i)); 
    //std::cout << " -- n: " << n << " -- gamma: " << gamma << n << " -- m(i): " << m(i)<< n << " -- K: " << K << " -- r(i): " << r(i) << " -- logr(i): " << logr(i) << std::endl;
  }

  Type likval;
  if (dbg > 0){
    std::cout << "" << std::endl;
    std::cout << "--- DEBUG: script start --- ans: " << ans << std::endl;
    for (int i=0; i < nstocks; i++){ 
      std::cout << "INPUT: logm(i): " << logm(i) << " -- i: " << i << std::endl; 
      std::cout << "INPUT: logn(i): " << logn(i) << std::endl;
      std::cout << "INPUT: logsdb(i): " << logsdb(i) << std::endl;
      std::cout << "INPUT: logK(i): " << logK(i) << std::endl;
    }
    for (int i=0; i < logphi.size(); i++){ 
      std::cout << "INPUT: logphi(i): " << logphi(i) << " -- i: " << i << std::endl; 
    }
    for (int i=0; i < logphipar.size(); i++){ 
      std::cout << "INPUT: logphipar(i): " << logphipar(i) << " -- i: " << i << std::endl; 
    }
    for (int i=0; i < nq; i++){ 
      std::cout << "INPUT: logq(i): " << logq(i) << " -- i: " << i << std::endl; 
    }
    for (int i=0; i < nqf; i++){ 
      std::cout << "INPUT: logsdf(i): " << logsdf(i) << " -- i: " << i << std::endl;
      std::cout << "INPUT: lambda(i): " << lambda(i) << std::endl;
    }
    std::cout << "INPUT: omega: " << omega << std::endl;

    std::cout << "logobsC.size(): " << logobsC.size() << "  Cpred.size(): " << Cpred.size() << "  logobsI.size(): " << logobsI.size() << "  dt.size(): " << dt.size() << "  logE.rows(): " << logE.rows() << "  logE.cols(): " << logE.cols() << "  logu.rows(): " << logu.rows() << "  logu.cols(): " << logu.cols() << "  B.rows(): " << B.rows() << "  B.cols(): " << B.cols() <<  "  mvec.size(): " << mvec.size() << "  iq.size(): " << iq.size() << "  ic.size(): " << ic.size() << "  logFmsy.size(): " << logFmsy.size() << "  logFmsyvec.size(): " << logFmsyvec.size() << "  logBmsy.size(): " << logBmsy.size() << "  logBmsyvec.size(): " << logBmsyvec.size() << "  m.size(): " << m.size() << "  logphi.size(): " << logphi.size() << "  logphipar.size(): " << logphipar.size() << std::endl;
    std::cout << "Bmsys: " << Bmsys << std::endl;
    std::cout << "Fmsys: " << Fmsys << std::endl;
    std::cout << "MSYs: " << MSYs << std::endl;
    std::cout << "Bmsyd: " << Bmsyd << std::endl;
    std::cout << "Fmsyd: " << Fmsyd << std::endl;
    std::cout << "MSYd: " << MSYd << std::endl;
  }
  // Calculate mvec if multiple rs are used (rarely the case).
  /*
  for (int i=0; i < ns; i++){
    ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
    if (dbg>1){
      std::cout << "-- i: " << i << " -- ind: " << ind << " -   mvec(i): " << mvec(i) << std::endl;
    }
    mvec(i) = m(ind);
    //logFmsyvec(i) = logFmsy(ind);
    //logBmsyvec(i) = logBmsy(ind);
  }
  */

  // PRIORS
  if (dbg > 0){
    std::cout << "PRIOR: priorn(0): " << priorn(0) << " -- priorn(1): " << priorn(1) << " -- priorn(2): " << priorn(2) << std::endl;
    std::cout << "PRIOR: priorr(0): " << priorr(0) << " -- priorr(1): " << priorr(1) << " -- priorr(2): " << priorr(2) << std::endl;
    std::cout << "PRIOR: priorK(0): " << priorK(0) << " -- priorK(1): " << priorK(1) << " -- priorK(2): " << priorK(2) << std::endl;
    std::cout << "PRIOR: priorm(0): " << priorm(0) << " -- priorm(1): " << priorm(1) << " -- priorm(2): " << priorm(2) << std::endl;
    std::cout << "PRIOR: priorq(0): " << priorq(0) << " -- priorq(1): " << priorq(1) << " -- priorq(2): " << priorq(2) << std::endl;
  }
  // Inverse gamma priors
  // TMB uses shape and scale parameterisation, but spict uses shape and rate similar to jags.
  // Prior for q. 
  if (prioriqgamma(2) == 1 & nq == 1){
    ans-= dgamma(1.0/exp(logq(0)), prioriqgamma(0), 1.0/prioriqgamma(1), 1); 
  }
  // Prior for sdb2. 
  if (priorisdb2gamma(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dgamma(1.0/sdb2(i), priorisdb2gamma(0), 1.0/priorisdb2gamma(1), 1); 
    }
  }
  // Prior for sdf2. 
  if (priorisdf2gamma(2) == 1){
    for (int i=0; i < nqf; i++){
      ans-= dgamma(1.0/sdf2(i), priorisdf2gamma(0), 1.0/priorisdf2gamma(1), 1); 
    }
  }
  // Prior for sdi2. 
  if (priorisdi2gamma(2) == 1){
    for (int i=0; i < nsdi; i++){
      ans-= dgamma(1.0/sdi2(i), priorisdi2gamma(0), 1.0/priorisdi2gamma(1), 1); 
    }
  }
  // Prior for sde2. 
  if (priorisde2gamma(2) == 1){
    for (int i=0; i < nsde; i++){
      ans-= dgamma(1.0/sde2(i), priorisde2gamma(0), 1.0/priorisde2gamma(1), 1); 
    }
  }
  // Log-normal priors
  if (priorn(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dnorm(logn(i), priorn(0), priorn(1), 1); // Prior for logn
    }
  }
  if (priorr(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dnorm(logr(i), priorr(0), priorr(1), 1); // Prior for logr
    }
  }
  if (priorK(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dnorm(logK(i), priorK(0), priorK(1), 1); // Prior for logK
    }
  }
  if (priorm(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dnorm(logm(i), priorm(0), priorm(1), 1); // Prior for logm
    }
  }
  if (priorq(2) == 1){
    for (int i=0; i < nq; i++){
      ans-= dnorm(logq(i), priorq(0), priorq(1), 1); // Prior for logq - log-normal
    }
  }
  if (priorqf(2) == 1){
    for (int i=0; i < nqf; i++){
      ans-= dnorm(logqf(i), priorqf(0), priorqf(1), 1); // Prior for logqf
    }
  }
  if (priorbkfrac(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dnorm(logB(i, 0) - logK(i), priorbkfrac(0), priorbkfrac(1), 1); // Prior for logbkfrac
    }
  }
  if (priorsdb(2) == 1){
    for (int i=0; i < nstocks; i++){
      ans-= dnorm(logsdb(i), priorsdb(0), priorsdb(1), 1); // Prior for logsdb
    }
  }
  if (priorsdf(2) == 1){
    for (int i=0; i < nqf; i++){
      ans-= dnorm(logsdf(i), priorsdf(0), priorsdf(1), 1); // Prior for logsdf
    }
  }
  if (priorsdi(2) == 1){
    for (int i=0; i < nsdi; i++){ 
      ans-= dnorm(logsdi(i), priorsdi(0), priorsdi(1), 1); // Prior for logsdi
    }
  }
  if (priorsde(2) == 1){
    for (int i=0; i < nsde; i++){ 
      ans-= dnorm(logsde(i), priorsde(0), priorsde(1), 1); // Prior for logsde
    }
  }
  if (priorsdc(2) == 1){
    for (int i=0; i < nqf; i++){ 
      ans-= dnorm(logsdc(i), priorsdc(0), priorsdc(1), 1); // Prior for logsdc
    }
  }
  if (prioralpha(2) == 1){
    for (int i=0; i < nsdi; i++){ 
      ans-= dnorm(logalpha(i), prioralpha(0), prioralpha(1), 1); // Prior for logalpha
    }
  }
  if (priorbeta(2) == 1){
    for (int i=0; i < nqf; i++){ 
      ans-= dnorm(logbeta(i), priorbeta(0), priorbeta(1), 1); // Prior for logbeta
    }
  }
  if (priorB(2) == 1 & nstocks == 1){
    ind = CppAD::Integer(priorB(4)-1);
    ans-= dnorm(logB(0, ind), priorB(0), priorB(1), 1); // Prior for logB
  }
  if (priorF(2) == 1 & nqf == 1){
    ind = CppAD::Integer(priorF(4)-1);
    ans-= dnorm(logF(0, ind), priorF(0), priorF(1), 1); // Prior for logF
  }
  if (priorBBmsy(2) == 1 & nstocks == 1){
    ind = CppAD::Integer(priorBBmsy(4)-1);
    ans-= dnorm(logB(0, ind) - logBmsyvec(ind), priorBBmsy(0), priorBBmsy(1), 1); // Prior for logBBmsy
  }
  if (priorFFmsy(2) == 1 & nqf == 1){
    ind = CppAD::Integer(priorFFmsy(4)-1);
    ans-= dnorm(logF(0, ind) - logFmsyvec(ind), priorFFmsy(0), priorFFmsy(1), 1); // Prior for logFFmsy
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
  //vector<Type> logFs = logF
  //vector<Type> logES(ns);
  matrix<Type> logES(nfleets, ns); // Seasonal component of E
  //matrix<Type> logE(nfleets, ns); // Seasonal component of E
  if (simple == 0){
    if (dbg > 0){
      std::cout << "--- DEBUG: F loop start --- ans: " << ans << std::endl;
    }
    for (int jj=0; jj < nfleets; jj++){
      // Diffusion component of E
      for (int i=1; i < ns; i++){
	Type Epredtmp;
	if (efforttype == 1.0){
	  Epredtmp = predictE1(logE(jj, i-1), dt(i), sdf2(jj));
	}
	if (efforttype == 2.0){
	  Epredtmp = predictE2(logE(jj, i-1), dt(i), sdf2(jj));
	}
	Type logEpred = log( ffacvec(i) * Epredtmp + fconvec(i) );;
	likval = dnorm(logE(jj, i), logEpred, sqrt(dt(i-1))*sdf(jj), 1);
	ans -= likval;
	// DEBUGGING
	if (dbg > 1){
	  std::cout << "-- i: " << i << "- jj: " << jj << " -   logE(jj, i-1): " << logE(jj, i-1) << "  logE(jj, i): " << logE(jj, i) << "  ffacvec(i): " << ffacvec(i) << "  fconvec(i): " << fconvec(i) << "  sdf(jj): " << sdf(jj) << "  efforttype: " << efforttype << "  likval: " << likval << "  ans:" << ans << std::endl;
	}
      }
      // Seasonal component
      if (dbg > 0){ 
	std::cout << "-- seasontype: " << seasontype << std::endl; 
      }
      for (int i=1; i < ns; i++){
	logES(jj, i) = 0.0; // Initialise
      }
      if (seasontype == 1.0){
	// Spline
	int ind2;
	for (int i=0; i < ns; i++){
	  ind2 = CppAD::Integer(seasonindex(i));
	  logES(jj, i) += seasonspline(ind2);
	  // DEBUGGING
	  if (dbg > 1){
	    //std::cout << "-- i: " << i << " -   logE(i): " << logE(i) << " logEs(i): " << logEs(i) << " ind2: " << ind2 << " seasonspline(ind2): " << seasonspline(ind2) << std::endl;
	    std::cout << "-- i: " << i << "- jj: " << jj << " -   logE(jj, i): " << logE(jj, i) << " logES(jj, i): " << logES(jj, i) << " ind2: " << ind2 << " seasonspline(ind2): " << seasonspline(ind2) << std::endl;
	  }
	}
      }
      if (seasontype == 2.0){
	// Coupled SDEs
	//for (int j=0; j < nsdu; j++){ // possible to use more than one constituent
	for (int j=0; j < 1; j++){ // Only using one constituent
	  Type per = j + 1.0;
	  if (dbg > 0){ 
	    std::cout << "-- j:" << j << "- per:" << per << "- omega:" << omega << std::endl; 
	  } 
	  for (int i=1; i < ns; i++){
	    // Analytical expression for matrix exponential
	    matrix<Type> trigmat(2, 2);
	    trigmat(0, 0) =  cos(per * omega * dt(i-1));
	    trigmat(0, 1) = -sin(per * omega * dt(i-1));
	    trigmat(1, 0) =  sin(per * omega * dt(i-1));
	    trigmat(1, 1) =  cos(per * omega * dt(i-1));
	    if (dbg > 0){ 
	      std::cout << "-- trigmat: " << trigmat << std::endl; 
	    }
	    matrix<Type> expmAt = exp(-lambda(jj) * dt(i-1)) * trigmat;
	    if (dbg > 0){ 
	      std::cout << "-- expmAt: " << expmAt << std::endl; 
	    }
	    // Corrected standard deviation
	    //Type sduana = sdu(j) * sqrt(1.0/(2.0*lambda(jj)) * (1.0 - exp(-2.0*lambda(jj)*dt(i-1))));
	    Type sduana = sdu(jj) * sqrt(1.0/(2.0*lambda(jj)) * (1.0 - exp(-2.0*lambda(jj)*dt(i-1))));
	    if (dbg > 0){ 
	      std::cout << "-- sduana: " << sduana << std::endl; 
	    }
	    // Extract master and slave states from logu
	    vector<Type> sublogumF = logu.col(i-1);
	    vector<Type> sublogum(2);
	    //sublogum(0) = sublogumF(2*j);   // Master state - old, for multiple constituents
	    //sublogum(1) = sublogumF(2*j+1); // Slave state - old, for multiple constituents
	    sublogum(0) = sublogumF(2*jj);   // Master state
	    sublogum(1) = sublogumF(2*jj+1); // Slave state
	    if (dbg > 0){ 
	      std::cout << "-- sublogumF: " << sublogumF << "- sublogum: " << sublogum << std::endl; 
	    }
	    vector<Type> logupred = expmAt * sublogum;
	    if (dbg > 0){ 
	      std::cout << "-- logupred: " << logupred << std::endl; 
	    }
	    likval = 0.0;
	    for (int k=0; k < logupred.size(); k++){ 
	      if (dbg > 0){ 
		std::cout << "-- k: " << k << "- 2*j+k: " << 2*j+k << " - logu(2*j+k, i): " << logu(2*j+k, i) << std::endl; 
	      }
	      //likval += dnorm(logu(2*j+k, i), logupred(k), sduana, 1); // old, for multiple constituents
	      likval += dnorm(logu(2*jj + k, i), logupred(k), sduana, 1); 
	    }
	    ans -= likval;
	    // DEBUGGING
	    if (dbg > 0){
	      std::cout << "-- i: " << i << "- jj: " << jj << " -   logu(0, i): " << logu(0, i) << "  logupred(0): " << logupred(0) << " -   logu(1, i): " << logu(1, i) << "  logupred(1): " << logupred(1) << "  sdu(jj): " << sdu(jj) << "  sduana: " << sduana << "  likval: " << likval << "  ans:" << ans << std::endl;
	    }
	  }
	  // Sum constituents
	  for (int i=0; i < ns; i++){
	    logES(jj, i) += logu(2*j, i); 
	  }
	}
      }
    }
  } else {
    for (int i=0; i < ns; i++){
      logES(i) = -30; // If using simple set fishing mortality to something small.
    }
  }
  // Sum diffusion and seasonal component
  matrix<Type> E(nfleets, ns);
  for (int jj=0; jj < nfleets; jj++){
    for (int i=0; i < ns; i++){
      E(jj, i) = exp(logES(jj, i) + logE(jj, i));
    }
  }


  //matrix<Type> logFs = log(F);

  // Sum Fs interacting with stock si
  matrix<Type> Fstock(nstocks, ns);
  for (int i=0; i < ns; i++){
    for (int si=0; si < nstocks; si++){
      Fstock(si, i) = 0.0; // Initialise
    }
    for (int k=0; k < nqf; k++){
      sind = CppAD::Integer(targetmap(k, 0) - 1); // minus 1 because R starts at 1 and c++ at 0
      //find = CppAD::Integer(targetmap(k, 1) - 1); // minus 1 because R starts at 1 and c++ at 0
      //Fstock(sind, i) += F(find, i);
      Fstock(sind, i) += F(k, i);
    }
  }

  // BIOMASS PREDICTIONS
  if (dbg > 0){
    std::cout << "--- DEBUG: B loop start --- ans: " << ans << std::endl;
  }
  matrix<Type> logBpred(nstocks, ns);
  for (int si=0; si < nstocks; si++){
    for (int i=0; i < (ns-1); i++){
      // To predict B(i) use dt(i-1), which is the time interval from t_i-1 to t_i
      if (simple == 0){
	logBpred(si, i+1) = predictlogB(B(si, i), Fstock(si, i), gamma(si), m(si), K(si), dt(i), n(si), sdb2(si));
      } /* else {
	Type Ftmp = 0.0;
	// Use naive approach
	Type Bpredtmp = exp(predictlogB(B(i), Ftmp, gamma(si), mvec(i), K, dt(i), n, sdb2) + 0.5*sdb2*dt(i)) - exp(logobsC(i));
	if (Bpredtmp < 0){
	  Bpredtmp = 1e-8; // Ugly ugly ugly hack to avoid taking log of negative
	}
	logBpred(i+1) = log(Bpredtmp);
	logFs(i) = logobsC(i) - logB(i); // Calculate fishing mortality
	} */
      likval = dnorm(logBpred(si, i+1), logB(si, i+1), sqrt(dt(i))*sdb(si), 1);
      ans -= likval;
      // DEBUGGING
      if (dbg > 1){
	std::cout << "-- i: " << i << "- si: " << si << " -   logB(si, i+1): " << logB(si, i+1) << "  logBpred(si, i+1): " << logBpred(si, i+1) << "  sdb(si): " << sdb(si) << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
    if (simple == 1){ 
      //logFs(ns-1) = logFs(ns-2);
    }
  }

  // EFFORT PREDICTIONS
  vector<Type> Ecumpred(nobsE);
  vector<Type> logEcumpred(nobsE);
  if (simple == 0){
    for (int i=0; i < nobsE; i++){
      Ecumpred(i) = 0.0;
      // Sum effort contributions from each sub interval
      for (int j=0; j < ne(i); j++){
	ind = CppAD::Integer(ie(i)-1) + j; // minus 1 because R starts at 1 and c++ at 0
	Ecumpred(i) += E(ind) * dt(ind);
      }
      logEcumpred(i) = log(Ecumpred(i));
      // DEBUGGING
      if (dbg>1){
	std::cout << "-- i: " << i << " -  ind: " << ind << "  logEcumpred(i): " << logEcumpred(i) << std::endl;
      }
    }
  }

  // CATCH PREDICTIONS
  matrix<Type> Cpredsub(nqf, ns);
  if (simple == 0){
    for (int k=0; k < nqf; k++){ 
      sind = CppAD::Integer(targetmap(k, 0) - 1); // minus 1 because R starts at 1 and c++ at 0
      find = CppAD::Integer(targetmap(k, 1) - 1); // minus 1 because R starts at 1 and c++ at 0
      for (int i=0; i < (ns-1); i++){ // ns-1 because dt is 1 shorter than state vec
	// For Cpredsub(i) use dt(i) because Cpredsub(i) is integrated over t_i to t_i+1
	Cpredsub(k, i) =  predictC(F(find, i), B(sind, i), dt(i));
      }
    }
    for (int i=0; i < nobsCp; i++){
      find = CppAD::Integer(iff(i) - 1); // minus 1 because R starts at 1 and c++ at 0
      // Sum catch contributions from each sub interval
      for (int j=0; j < nc(i); j++){
	ind = CppAD::Integer(ic(i) - 1) + j; // minus 1 because R starts at 1 and c++ at 0
	Cpred(i) += Cpredsub(find, ind);
      }
      logCpred(i) = log(Cpred(i));
      // DEBUGGING
      if (dbg > 1){
	std::cout << "-- i: " << i << " -  ind: " << ind << "  logCpred(i): " << logCpred(i) << std::endl;
      }
    }
  } else {
    for (int i=0; i < (ns-1); i++){ // ns-1 because logobsC is 1 shorter than ns
      //Cpredsub(i) = exp(logobsC(i));
      //logCpred(i) = logobsC(i);
    }
  }

  // CALCULATE PRODUCTION
  matrix<Type> Cpredsubperstock(nstocks, ns);
  for (int i=0; i < (ns-1); i++){
    for (int si=0; si < nstocks; si++){
      Cpredsubperstock(si, i) = 0.0; // Initialise
    }
    for (int k=0; k < nqf; k++){
      sind = CppAD::Integer(targetmap(k, 0) - 1); // minus 1 because R starts at 1 and c++ at 0
      Cpredsubperstock(sind, i) += Cpredsub(k, i);
    }
  }
  matrix<Type> P(nstocks, ns-1);
  for (int si=0; si < nstocks; si++){
    for (int i=0; i < (ns-1); i++){
      P(si, i) = B(si, i+1) - B(si, i) + Cpredsubperstock(si, i);
    }
  }

  /*
    --- OBSERVATION EQUATIONS ---
  */
  int inds;

  // CATCHES
  int indsdc;
  int indscpred;
  if (simple == 0){
    if (dbg > 0){
      std::cout << "--- DEBUG: Cpred loop start --- ans: " << ans << std::endl;
    }
    // fac and pp are used for the outlier robust Gaussian mixture.
    for (int i=0; i < nobsC; i++){
      inds = CppAD::Integer(isc(i)-1);
      indsdc = CppAD::Integer(isdc(i)-1);
      indscpred = CppAD::Integer(icpred(i)-1);
      if (robflagc == 1.0){
	likval = log(pp*dnorm(logCpred(indscpred), logobsC(i), stdevfacc(i)*sdc(indsdc), 0) + (1.0-pp)*dnorm(logCpred(indscpred), logobsC(i), robfac*stdevfacc(i)*sdc(indsdc), 0));
      } else {
	likval = dnorm(logCpred(indscpred), logobsC(i), stdevfacc(i)*sdc(indsdc), 1);
      }
      ans-= keep(inds) * likval;
      // DEBUGGING
      if (dbg > 1){
	std::cout << "-- i: " << i << "- indscpred: " << indscpred << " -   logobsC(i): " << logobsC(i) << " -   stdevfacc(i): " << stdevfacc(i) << "  sdc(indsdc): " << sdc(indsdc) << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
  }

  // EFFORT
  int indsde;
  if (simple == 0){
    if (dbg > 0){
      std::cout << "--- DEBUG: Epred loop start --- ans: " << ans << std::endl;
    }
    for (int i=0; i < nobsE; i++){
      indsde = CppAD::Integer(isde(i)-1);
      //logEpred(i) = logFcumpred(i) - logqf; // E = 1/q * integral{F_t dt}
      logEpred(i) = logEcumpred(i); // E = 1/q * integral{F_t dt}
      if (robflage == 1.0){
	likval = log(pp*dnorm(logEpred(i), logobsE(i), stdevface(i)*sde(indsde), 0) + (1.0-pp)*dnorm(logEpred(i), logobsE(i), robfac*stdevface(i)*sde(indsde), 0));
      } else {
	likval = dnorm(logEpred(i), logobsE(i), stdevface(i)*sde(indsde), 1);
      }
      inds = CppAD::Integer(ise(i)-1);
      ans-= keep(inds) * likval;
      // DEBUGGING
      if (dbg > 1){
	std::cout << "-- i: " << i << " -   logobsE(i): " << logobsE(i) << " -   stdevface(i): " << stdevface(i) << "  sde(indsde): " << sde(indsde) << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
  }

  // BIOMASS INDEX
  if (dbg > 0){
    std::cout << "--- DEBUG: Ipred loop start --- ans: " << ans << std::endl;
    std::cout << " logobsI.size(): " << logobsI.size() << "  iq.size(): " << iq.size() << "  ii.size(): " << ii.size() << "  logq.size(): " << logq.size() << "  nobsI: " << nobsI << std::endl;
  }
  int indb;
  int indq;
  int indsdi;
  for (int i=0; i < nobsI; i++){
    ind = CppAD::Integer(ii(i)-1);
    indb = CppAD::Integer(ib(i)-1);
    indq = CppAD::Integer(iq(i)-1);
    indsdi = CppAD::Integer(isdi(i)-1);
    inds = CppAD::Integer(isi(i)-1);
    logIpred(i) = logq(indq) + log(B(indb, ind));
    if (robflagi == 1.0){
      likval = log(pp*dnorm(logobsI(i), logIpred(i), stdevfaci(i)*sdi(indsdi), 0) + (1.0-pp)*dnorm(logobsI(i), logIpred(i), robfac*stdevfaci(i)*sdi(indsdi), 0));
    } else {
      likval = dnorm(logobsI(i), logIpred(i), stdevfaci(i)*sdi(indsdi), 1);
    }
    ans -= keep(inds) * likval;
    // DEBUGGING
    if (dbg > 1){
      std::cout << "-- i: " << i << " -  ind: " << ind << " -  indq: " << indq << " -  indb: " << indb << " -  indsdi: " << indsdi << " -  inds: " << inds << " -   logobsI(i): " << logobsI(i) << "  logIpred(i): " << logIpred(i) << "  stdevfaci(i): " << stdevfaci(i) << "  likval: " << likval << "  sdi: " << sdi << "  ans:" << ans << std::endl;
    }
  }



  /*
  --- ONE-STEP-AHEAD PREDICTIONS ---
  */

  if (dbg > 0){
    std::cout << "--- DEBUG: ONE-STEP-AHEAD PREDICTIONS --- ans: " << ans << std::endl;
    std::cout << "-- dtpredcnsteps: " << dtpredcnsteps << "  dtpredcinds.size(): " << dtpredcinds.size() <<std::endl;
    std::cout << "-- dtpredensteps: " << dtpredensteps << "  dtpredeinds.size(): " << dtpredeinds.size() <<std::endl;
  }
  // Calculate logFstock
  matrix<Type> logFstock(nstocks, ns);
  for (int i=0; i < ns; i++){
    for (int si=0; si < nstocks; si++){
      logFstock(si, i) = log(Fstock(si, i));
    }
  }

  // Catch predictions
  vector<Type> Cp(nqf);
  vector<Type> Cstockp(nstocks);
  for (int si=0; si < nstocks; si++){
    Cstockp(si) = 0.0;
  }
  for (int k=0; k < nqf; k++){
    Cp(k) = 0.0;
    for (int i=0; i < dtpredcnsteps; i++){
      ind = CppAD::Integer(dtpredcinds(i) - 1);
      if (dbg > 1){
	std::cout << "-- i: " << i << "-- k: " << k << " -  dtpredcinds(i)-1: " << ind << std::endl;
      }
      Cp(k) += Cpredsub(k, ind);
    }
    sind = CppAD::Integer(targetmap(k, 0) - 1);
    Cstockp(sind) += Cp(k);
  }
  vector<Type> logCp = log(Cp);
  vector<Type> logCstockp = log(Cstockp);

  // Effort prediction
  vector<Type> Ep(nfleets);
  for (int k=0; k < nfleets; k++){
    for (int i=0; i < dtpredensteps; i++){
      ind = CppAD::Integer(dtpredeinds(i) - 1);
      if (dbg > 1){
	std::cout << "-- i: " << i << " -  dtpredeinds(i)-1: " << ind << std::endl;
      }
      Ep(k) += E(k, ind);
    }
  }
  vector<Type> logEp = log(Ep);

  // Biomass and F at the end of the prediction time interval
  int pind = CppAD::Integer(dtprediind - 1);
  vector<Type> Bp = B.col(pind); 
  vector<Type> logBp = log(Bp);
  vector<Type> logBpBmsy(nstocks);
  vector<Type> logFstockp = logFstock.col(pind);
  vector<Type> logFstockpFmsy(nstocks);
  for (int si=0; si < nstocks; si++){ 
    logBpBmsy(si) = logBp(si) - logBmsyvec(si, pind);
    logFstockpFmsy(si) = logFstockp(si) - logFmsyvec(si, pind);
  }
  vector<Type> logBpK = logBp - logK;
  vector<Type> logFp = logF.col(pind); 
  vector<Type> logFpFmsy(nqf);
  for (int k=0; k < nqf; k++){
    sind = CppAD::Integer(targetmap(k, 0) - 1);
    logFpFmsy(k) = logFp(k) - logFmsyvec(sind, pind);
  }
  vector<Type> logIp(nindex);
  for (int i=0; i < nindex; i++){
    int sind = CppAD::Integer(index2sdb(i) - 1);
    int qind = CppAD::Integer(index2q(i) - 1);
    logIp(i) = logq(qind) + log(Bp(sind));
  }

  // Biomass and fishing mortality at last time point
  int lind = indlastobs - 1;
  vector<Type> logBl = logB.col(lind);
  vector<Type> logBlBmsy(nstocks);
  vector<Type> logFstockl = logFstock.col(lind);
  vector<Type> logFstocklFmsy(nstocks);
  for (int si=0; si < nstocks; si++){ 
    logBlBmsy(si) = logBl(si) - logBmsyvec(si, lind);
    logFstocklFmsy(si) = logFstockl(si) - logFmsyvec(si, lind);
  }
  vector<Type> logBlK = logBl - logK;
  vector<Type> logFl = logF.col(lind);
  vector<Type> logFlFmsy(nqf);
  for (int k=0; k < nqf; k++){
    sind = CppAD::Integer(targetmap(k, 0) - 1);
    logFlFmsy(k) = logFl(k) - logFmsyvec(sind, lind);
  }

  // Calculate relative levels of biomass and fishing mortality per stock
  matrix<Type> logBBmsy(nstocks, ns);
  matrix<Type> logFstockFmsy(nstocks, ns);
  for (int si=0; si < nstocks; si++){ 
    for (int i=0; i < ns; i++){ 
      logBBmsy(si, i) = logB(si, i) - logBmsyvec(si, i); 
      logFstockFmsy(si, i) = logFstock(si, i) - logFmsyvec(si, i); 
    }
  }
  // Calculate relative fishing mortality per fishery
  matrix<Type> logFFmsy(nqf, ns);
  for (int k=0; k < nqf; k++){ 
    sind = CppAD::Integer(targetmap(k, 0) - 1); // minus 1 because R starts at 1 and c++ at 0
    for (int i=0; i < ns; i++){ 
      logFFmsy(k, i) = logF(k, i) - logFmsyvec(sind, i); 
    }
  }

  // Calculate logbkfrac i.e. log[B(0)/K]
  vector<Type> logbkfrac(nstocks);
  for (int si=0; si < nstocks; si++){ 
    logbkfrac(si) = logB(si, 0) - logK(si);
  }

  // ADREPORTS
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
  ADREPORT(logFstockp);
  ADREPORT(logFstockpFmsy);
  ADREPORT(logFl);
  ADREPORT(logFlFmsy);
  ADREPORT(logFstockl);
  ADREPORT(logFstocklFmsy);
  ADREPORT(MSY);
  ADREPORT(MSYd);
  ADREPORT(MSYs);
  ADREPORT(logMSY);
  ADREPORT(logMSYd);
  ADREPORT(logMSYs);
  /*
  ADREPORT(Emsy);
  ADREPORT(Emsy2);
  ADREPORT(logEmsy);
  ADREPORT(logEmsy2);
  */
  ADREPORT(logbkfrac);
  ADREPORT(seasonsplinefine);
  // PREDICTIONS
  ADREPORT(Cp);
  ADREPORT(logIp);
  ADREPORT(logCp);
  ADREPORT(logCstockp);
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
  if (reportall){ 
    // These reports are derived from the random effects and are therefore vectors. TMB calculates the covariance of all sdreports leading to a very large covariance matrix which may cause memory problems.
    // B
    ADREPORT(logBBmsy);
    // F
    ADREPORT(logF);
    ADREPORT(logFstock);
    ADREPORT(logFFmsy);
    ADREPORT(logFstockFmsy);
    // C
    ADREPORT(logCpred);
    // I
    ADREPORT(logIpred);
    // E
    ADREPORT(logEpred);
  }

  // REPORTS (these don't require sdreport to be output)
  REPORT(Cp);
  REPORT(Ep);
  REPORT(logIp);
  REPORT(MSY);
  REPORT(Bmsy);
  REPORT(Fmsy);
  REPORT(stochmsy);
  REPORT(P);

  return ans;
}

