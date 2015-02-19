// Surplus Production in Continuous-Time (SPiCT)
// 14.10.2014
#include <TMB.hpp>

/* t-distribution */
template<class Type>
Type ltdistr(const Type &x, const Type &df)
{
  Type p = 1.0;
  Type h = 0.5;
  return -(-lgamma(h*(df+p)) + h*log(df*M_PI) + lgamma(h*df) + h*(df+p)*log(Type(1.0)+x*x/df));
}

/* Calculate Binf */
template<class Type>
Type calculateBinf(const Type &K, const Type &F, const Type &gamma, const Type &m, const Type &n, const Type &sdb2)
{
  return K*pow(1.0 - K/(gamma*m)*F - K/(gamma*m)*0.5*sdb2, 1.0/(n-1.0));
}

/* Predict log F */
template<class Type>
Type predictlogF(const Type &phi1, const Type &logF1, const Type &phi2, const Type &logF2)
{
  return phi1*logF1 + phi2*logF2;
}

/* Predict biomass */
template<class Type>
Type predictB(const Type &B0, const Type &F, const Type &gamma, const Type &m, const Type &K, const Type &dt, const Type &n, const Type &sdb2=0)
{
  // Pella-Tomlinson
  return exp( log(B0) + (gamma*m/K - gamma*m/K*pow(B0/K, n-1.0) - F - 0.5*sdb2)*dt ); // Euler
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

  DATA_INTEGER(reportall); // Euler flag.
  DATA_VECTOR(dt);         // Time steps
  DATA_VECTOR(dtpredcinds);
  DATA_INTEGER(dtpredcnsteps); // Number of sub time step for prediction
  DATA_SCALAR(dtprediind);
  DATA_INTEGER(indlastobs);// Index of B and F corresponding to the last observation.
  DATA_VECTOR(obsC);       // Catches
  DATA_VECTOR(ic);         // Vector such that B(ii(i)) is the state corresponding to obsC(i)
  DATA_VECTOR(nc);         // nc(i) gives the number of time intervals obsC(i) spans
  DATA_VECTOR(I);          // Index
  DATA_VECTOR(ii);         // A vector such that B(ii(i)) is the state corresponding to I(i)
  DATA_VECTOR(iq);         // A vector such that iq(i) is the index number corresponding to I_iq(i)
  DATA_VECTOR(ir);         // A vector indicating when the different rs should be used
  DATA_VECTOR(seasons);    // A vector of length ns indicating to which season a state belongs
  DATA_VECTOR(seasonindex);// A vector of length ns giving the number stepped within the current year
  DATA_MATRIX(splinemat);  // Design matrix for the seasonal spline
  DATA_SCALAR(ffac);       // Management factor each year multiply the predicted F with ffac
  DATA_VECTOR(indpred);    // A vector indicating when the management factor should be applied
  DATA_SCALAR(robflagc);        // Catch Degrees of freedom of t-distribution (only used if tdf < 25)
  DATA_SCALAR(robflagi);        // Index Degrees of freedom of t-distribution (only used if tdf < 25)
  DATA_INTEGER(lamperti);  // Lamperti flag.
  DATA_INTEGER(euler);     // Euler flag.
  DATA_INTEGER(stochmsy);  // Use stochastic msy?
  DATA_SCALAR(dbg);        // Debug flag, if == 1 then print stuff.
  PARAMETER_VECTOR(logphi);   // Season levels of F.
  PARAMETER(logitpp);      // 
  PARAMETER(logp1robfac);  // 
  PARAMETER(logalpha);     // sdi = alpha*sdb
  PARAMETER(logbeta);      // sdc = beta*sdf
  PARAMETER(logbkfrac);    // B0/K fraction
  //PARAMETER(logF0);        // F at time 0
  PARAMETER_VECTOR(logm);  // m following the Fletcher formulation (see Prager 2002)
  PARAMETER(logK);         // Carrying capacity
  PARAMETER_VECTOR(logq);  // Catchability
  PARAMETER(logn);         // Pella-Tomlinson exponent
  PARAMETER(logsdf);       // Standard deviation for F
  PARAMETER(logsdb);       // Standard deviation for Index
  PARAMETER_VECTOR(logF);  // Random effects vector
  PARAMETER_VECTOR(logB);  // Random effects vector
  PARAMETER(dum);       // Dummy parameter, needed because the RE cannot be evaluated without FE
  ans+=dum*dum; // Add dummy contribution (0 for dum=0)

  lamperti = 1.0; // Not used anymore
  euler = 1.0; // Not used anymore
  vector<Type> logphipar(logphi.size()+1);
  logphipar(0) = 0.0; // The first logphi is set to 0, the rest are estimated relative to this.
  for(int i=1; i<logphipar.size(); i++){ logphipar(i) = logphi(i-1); }
  vector<Type> temp = splinemat.col(0);
  vector<Type> seasonspline(temp.size());
  seasonspline = splinemat * logphipar;
  Type pp = 1.0/(1.0 + exp(-logitpp));
  Type robfac = 1.0 + exp(logp1robfac);
  Type bkfrac = exp(logbkfrac);
  //Type F0 = exp(logF0);
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
  Type sdf = exp(logsdf);
  Type sdb = exp(logsdb);
  Type sdb2 = sdb*sdb;
  Type sdi = exp(logalpha)*sdb;
  Type sdc = exp(logbeta)*sdf;
  Type logsdi = log(sdi);
  Type logsdc = log(sdc);
  int nobsC = obsC.size();
  int nobsCp = ic.size();
  int nobsI = I.size();
  int ns = logF.size();
  vector<Type> P(ns-1);
  vector<Type> Binf(ns);
  vector<Type> logBinf(ns);
  vector<Type> B = exp(logB);
  Type logB0 = logbkfrac + logK;
  vector<Type> Bpred(ns);
  vector<Type> mvec(ns);
  vector<Type> ffacvec(ns);
  for(int i=0; i<ns; i++){ ffacvec(i) = 1.0; }
  vector<Type> Cpred(nobsCp);
  for(int i=0; i<nobsCp; i++){ Cpred(i) = 0.0; }
  vector<Type> Cpredsub(ns);
  vector<Type> logIpred(nobsI);
  vector<Type> logCpred(nobsCp);
  Type mmean = m(0);
  if(nm > 1) for(int i=1; i<nm; i++){ mmean += m(i)/nm; }

  // Deterministic reference points
  Type Bmsyd = K * pow(1.0/n, 1.0/(n-1.0));
  //Type MSYd = mmean;
  Type MSYd = m(0);
  Type Fmsyd = MSYd/Bmsyd;
  Type logBmsyd = log(Bmsyd);
  Type logMSYd = log(MSYd);
  Type logFmsyd = log(Fmsyd);

  // Stochastic reference points (NOTE: only proved for n>1, Bordet and Rivest (2014))
  Type rbmean = (n-1)/n*gamma*mmean/K;
  Type Bmsys = Bmsyd * (1.0 - (1.0 + rbmean*(p-1.0)/2.0)*sdb2 / (rbmean*pow(2.0-rbmean, 2.0)));
  if(dbg==-1){
    std::cout << "--- DEBUG Bmsy ---" << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
    std::cout << "mmean: " << mmean << std::endl;
    std::cout << "K: " << K << std::endl;
    std::cout << "Bmsyd: " << Bmsyd << std::endl;
    std::cout << "rbmean: " << rbmean << std::endl;
    std::cout << "sdb2: " << sdb2 << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "Bmsys: " << Bmsys << std::endl;
  }
  Type Fmsys = Fmsyd - (p*(1.0-rbmean)*sdb2) / pow(2.0-rbmean, 2.0) ;
  Type MSYs = MSYd * (1.0 - ((p+1.0)/2.0*sdb2) / (1.0 - pow(1.0-rbmean, 2.0)));
  Type logBmsys = log(Bmsys);
  Type logFmsys = log(Fmsys);
  Type logMSYs = log(MSYs);

  Type Bmsy;
  Type MSY;
  Type Fmsy;
  Type logBmsy;
  Type logFmsy;
  Type logMSY;

  if(stochmsy==1){
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
    Emsy(i) = Fmsy/q(i); 
    Emsy2(i) = Fmsy/exp(logq2(i))*1.0e4; 
  }

  // Calculate growth rate
  Type sign = 1.0;
  if(n < 1) sign = -1.0; // Following Fletcher (1978)
  vector<Type> r(nm);
  vector<Type> logr(nm);
  for(int i=0; i<nm; i++){ 
    r(i) = sign * gamma * m(i) / K; 
    logr(i) = log(r(i)); 
  }

  Type likval;

  if(dbg>0){
    std::cout << "--- DEBUG: script start ---" << std::endl;
    //for(int i=0; i<ns; i++) std::cout << "F(i): " << F(i) << std::endl;
    //std::cout << "INPUT: logbkfrac: " << logbkfrac << std::endl;
    //std::cout << "INPUT: logF0: " << logF0 << std::endl;
    for(int i=0; i<nm; i++){ std::cout << "INPUT: logm(i): " << logm(i) << " -- i: " << i << std::endl; }
    for(int i=0; i<logphi.size(); i++){ std::cout << "INPUT: logphi(i): " << logphi(i) << " -- i: " << i << std::endl; }
    for(int i=0; i<logphipar.size(); i++){ std::cout << "INPUT: logphipar(i): " << logphipar(i) << " -- i: " << i << std::endl; }
    std::cout << "INPUT: logK: " << logK << std::endl;
    for(int i=0; i<nq; i++){ std::cout << "INPUT: logq(i): " << logq(i) << " -- i: " << i << std::endl; }
    std::cout << "INPUT: logn: " << logn << std::endl;
    std::cout << "INPUT: logsdf: " << logsdf << std::endl;
    std::cout << "INPUT: logsdb: " << logsdb << std::endl;
    std::cout << "obsC.size(): " << obsC.size() << "  Cpred.size(): " << Cpred.size() << "  I.size(): " << I.size() << "  dt.size(): " << dt.size() << "  logF.size(): " << logF.size() << "  B.size(): " << B.size() << "  P.size(): " << P.size() << "  mvec.size(): " << mvec.size() << "  iq.size(): " << iq.size() << "  ic.size(): " << ic.size() << "  logphi.size(): " << logphi.size() << "  logphipar.size(): " << logphipar.size() << std::endl;
  }
  // Calculate mvec if multiple rs are used (rarely the case).
  int ind;
  for(int i=0; i<ns; i++){
    ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
    mvec(i) = m(ind);
    if(dbg>1){
      std::cout << "-- i: " << i << "-- ind: " << ind << " -   mvec(i): " << mvec(i) << std::endl;
    }
  }
  for(int i=1; i<indpred.size(); i++){ // don't use i=0 because this is only for plotting
    ind = CppAD::Integer(indpred(i)-1); // minus 1 because R starts at 1 and c++ at 0
    ffacvec(ind) = ffac;
    if(dbg>1){
      std::cout << "-- i: " << i << "-- ind: " << ind << " -   ffacvec(i): " << ffacvec(i) << std::endl;
    }
  }
  //  for(int i=0; i<ns; i++){
  //  if(F(i)==mvec(i)) std::cout << "Warning: F(i)-mvec(i): " << F(i)-mvec(i) << std::endl;
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
  Type sd = 1e-6; // This one is used in the hack to fix B0 and F0.

  // FISHING MORTALITY
  if(dbg>0){
    std::cout << "--- DEBUG: F loop start" << std::endl;
  }
  // Hack to set log(F(0)) equal to the fixed effect log(F0).
  //likval = dnorm(logF0, logF(0), sd, 1);
  //ans-=likval;
  if(dbg>0){
    //std::cout << "- logF0: " << logF0 << " logF(0): " << logF(0) << "  likval: " << likval << "  ans:" << ans <<  std::endl;
  }
  for(int i=1; i<ns; i++){
    //Type logFpred = log(ffacvec(i)) + predictlogF(phi1, logF(i-1), phi2, logF(i-delay));
    Type logFpred = log(ffacvec(i)) + logF(i-1);
    likval = dnorm(logF(i), logFpred, sqrt(dt(i-1))*sdf, 1);
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logF(i-1): " << logF(i-1) << "  logF(i): " << logF(i) << "  ffacvec(i): " << ffacvec(i) << "  sdf: " << sdf << "  likval: " << likval << "  ans:" << ans << std::endl;
    }
  }
  vector<Type> logFs(ns);

  int ind2;
  for(int i=0; i<ns; i++){
    ind2 = CppAD::Integer(seasonindex(i));
    logFs(i) = seasonspline(ind2) + logF(i);
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logF(i): " << logF(i) << " logFs(i): " << logFs(i) << " ind2: " << ind2 << " seasonspline(ind2): " << seasonspline(ind2) << std::endl;
    }
  }
  vector<Type> F = exp(logFs);

  // CALCULATE B_infinity
  for(int i=0; i<ns; i++) Binf(i) = K * pow( 1.0 - (n-1.0)*F(i)/(n * Fmsyd), 1.0/(n-1.0)) * (1 - n/2.0 / (1.0- pow(1.0-n*Fmsyd + (n-1.0)*F(i), 2.0))* sdb2); // Bordet & Rivest (2014) eqn. 9
  logBinf = log(Binf);

  // BIOMASS PREDICTIONS
  // Hack to set log(B(0)) equal to the fixed effect log(B0).
  likval = dnorm(logB0, log(B(0)), sd, 1);
  ans-=likval;
  if(dbg>1){
    std::cout << "-- i: " << 0 << " -   logB0: " << logB0 << "  log(B(0)): " << log(B(0)) << "  sd: " << sd << "  likval: " << likval << "  ans:" << ans << std::endl;
  }
  for(int i=0; i<(ns-1); i++){
    // To predict B(i) use dt(i-1), which is the time interval from t_i-1 to t_i
    Bpred(i+1) = predictB(B(i), F(i), gamma, mvec(i), K, dt(i), n, sdb2);
    likval = dnorm(log(Bpred(i+1)), logB(i+1), sqrt(dt(i))*sdb, 1);
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logB(i+1): " << logB(i+1) << "  log(Bpred(i+1)): " << log(Bpred(i+1)) << "  sdb: " << sdb << "  likval: " << likval << "  ans:" << ans << std::endl;
    }
  }

  // CATCH PREDICTIONS
  for(int i=0; i<(ns-1); i++){ // ns-1 because dt is 1 shorter than state vec
    // For Cpredsub(i) use dt(i) because Cpredsub(i) is integrated over t_i to t_i+1
    Cpredsub(i) =  predictC(F(i), B(i), dt(i));
  }
  for(int i=0; i<nobsCp; i++){
    // Sum catch contributions from each sub interval
    for(int j=0; j<nc(i); j++){
      ind = CppAD::Integer(ic(i)-1) + j; // minus 1 because R starts at 1 and c++ at 0
      Cpred(i) += Cpredsub(ind);
      logCpred(i) = log(Cpred(i));
    }
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -  ind: " << ind << "  logCpred(i): " << logCpred(i) << std::endl;
    }
  }

  // CALCULATE PRODUCTION
  for(int i=0; i<(ns-1); i++) P(i) = B(i+1) - B(i) + Cpredsub(i);



  /*
  --- OBSERVATION EQUATIONS ---
  */

  // CATCHES
  if(dbg>0){
    std::cout << "--- DEBUG: Cpred loop start" << std::endl;
  }
  // fac and pp are used for the outlier robust Gaussian mixture.
  //Type robfac = 20.0;
  //Type pp = 0.95;
  for(int i=0; i<nobsC; i++){
    int j = CppAD::Integer(nc(i)-1);
    ind = CppAD::Integer(ic(i)-1) + j; // minus 1 because R starts at 1 and c++ at 0
    if(robflagc==1.0){
      //Type z = (logCpred(i)-log(obsC(i)))/sdc;
      //likval = -log(sdc) + dnorm(z, Type(0.0), Type(1.0), 1);
      //likval = -log(sdc) + ltdistr(z, Type(100.0));
      likval = log(pp*dnorm(logCpred(i), log(obsC(i)), sdc, 0) + (1.0-pp)*dnorm(logCpred(i), log(obsC(i)), robfac*sdc, 0));
    } else {
      likval = dnorm(logCpred(i), log(obsC(i)), sdc, 1);
    }
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logobsC(i): " << log(obsC(i))<< "  sdc: " << sdc << "  likval: " << likval << "  ans:" << ans << std::endl;
    }
  }

  // BIOMASS INDEX
  if(dbg>0){
    std::cout << "--- DEBUG: Ipred loop start" << std::endl;
    std::cout << " I.size(): " << I.size() << "  iq.size(): " << iq.size() << "  ii.size(): " << ii.size() << "  logq.size(): " << logq.size() << std::endl;
  }
  int indq;
  for(int i=0; i<nobsI; i++){
    if(I(i)>0){
      ind = CppAD::Integer(ii(i)-1);
      indq = CppAD::Integer(iq(i)-1);
      logIpred(i) = logq(indq) + log(B(ind));
      if(robflagi==1.0){
	//Type z = (log(I(i)) - logIpred(i))/sdi;
	//likval = -log(sdi) + dnorm(z, Type(0.0), Type(1.0), 1);
	//likval = -log(sdi) + ltdistr(z, Type(100.0));
	likval = log(pp*dnorm(log(I(i)), logIpred(i), sdi, 0) + (1.0-pp)*dnorm(log(I(i)), logIpred(i), robfac*sdi, 0));
      } else {
	likval = dnorm(log(I(i)), logIpred(i), sdi, 1);
      }
      ans-=likval;
      // DEBUGGING
      if(dbg>1){
	std::cout << "-- i: " << i << " -  ind: " << ind << " -  indq: " << indq << " -   log(I(i)): " << log(I(i)) << "  logIpred(i): " << logIpred(i) << "  sdi: " << sdi << "  likval: " << likval << "  ans:" << ans << std::endl;
      }
    }
  }



  /*
  --- ONE-STEP-AHEAD PREDICTIONS ---
  */

  if(dbg>0){
    std::cout << "--- DEBUG: ONE-STEP-AHEAD PREDICTIONS" << std::endl;
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
  Type Bp = B(CppAD::Integer(dtprediind-1)); 
  Type logBp = log(Bp);
  Type logBpBmsy = logBp - logBmsyd;
  Type logBpK = logBp - logK;
  Type logFp = logFs(CppAD::Integer(dtprediind-1)); 
  Type logFpFmsy = logFp - logFmsyd;

  vector<Type> logIp(nq);
  for(int i=0; i<nq; i++){
    logIp(i) = logq(i) + log(Bp);
  }

  // Biomass and fishing mortality at last time point
  // dtpredcinds(0) is the index of B and F corresponding to the time of the last observation.
  Type logBl = logB(indlastobs-1);
  Type logBlBmsy = logBl - logBmsy;
  Type logBlK = logBl - logK;
  Type logFl = logFs(indlastobs-1);
  Type logFlFmsy = logFl - logFmsy;

  // Biomass and fishing mortality over msy levels
  vector<Type> logBBmsy(ns);
  vector<Type> logFFmsy(ns);
  for(int i=0; i<ns; i++){ 
    logBBmsy(i) = logB(i) - logBmsy; 
    logFFmsy(i) = logFs(i) - logFmsy; 
  }

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
  // PREDICTIONS
  ADREPORT(Cp);
  ADREPORT(logIp);
  ADREPORT(logCp);
  // PARAMETERS
  ADREPORT(r);
  ADREPORT(logr);
  ADREPORT(K);
  ADREPORT(q);
  ADREPORT(logq2);
  ADREPORT(p);
  ADREPORT(gamma);
  ADREPORT(m);
  ADREPORT(sdf);
  ADREPORT(sdc);
  ADREPORT(sdb);
  ADREPORT(sdi);
  ADREPORT(logsdi);
  ADREPORT(logsdc);
  if(reportall){ 
    // These reports are derived from the random effects and are therefore vectors. TMB calculates the covariance of all sdreports leading to a very large covariance matrix which may cause memory problems.
    // B
    ADREPORT(logBinf);
    ADREPORT(logBBmsy);
    vector<Type> Bind(ns);
    for(int i=0; i<ns; i++) Bind(i) = 1.0 - B(i)/Binf(i);
    ADREPORT(Bind);
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

  return ans;
}

