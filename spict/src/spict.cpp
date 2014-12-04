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

/* Predict log F */
template<class Type>
Type predictlogF(const Type &phi1, const Type &logF1, const Type &phi2, const Type &logF2)
{
  return phi1*logF1 + phi2*logF2;
}

/* Calculate B_infinity */
template<class Type>
Type calculateBinf(const Type &K, const Type &F, const Type &r, const Type &p, const Type &sdb2=0, int lamperti=0)
{
  if(lamperti){
    // From Bordet & Rivest (2014)
    return K * pow(1.0 - F/r, 1.0/p) * (1.0 - (p+1.0)/2.0 / (1.0-pow(1.0-p*r+p*F, 2.0)) * sdb2);
  } else {
    return K * pow(1.0 - F/r, 1.0/p);
  }
}

/* Predict biomass */
template<class Type>
Type predictB(const Type &B0, const Type &Binf, const Type &F, const Type &r, const Type &K, const Type &dt, const Type &p, const Type &sdb2=0, int lamperti=0, int euler=0)
{
  if(euler) lamperti = 1;
  Type rate;
  if(lamperti){
    rate = r - F - 0.5*sdb2;
  } else {
    rate = r - F;    
  }
  if(euler){
    // Pella-Tomlinson
    return exp( log(B0) + (rate - r*pow(B0/K, p))*dt ); // Euler
  } else {
    // Schaefer only
    return 1 / ( 1/Binf + (1/B0 - 1/Binf) * exp(-rate*dt) ); // Approximative analytical, p=1
  }
}

/* Predict  catch*/
template<class Type>
Type predictC(const Type &F, const Type &K, const Type &r, const Type &B0, const Type &Binf, const Type &dt, const Type &sdb2=0, int lamperti=0, int euler=0)
{
  if(euler) lamperti = 1;
  Type rate;
  if(lamperti){
    rate = r - F - 0.5*sdb2;
  } else {
    rate = r - F;    
  }
  if(euler){
    return F*B0*dt;
  } else {
    return K/r*F * log( 1 - B0/Binf * (1 - exp(rate*dt)));
  }
}

/* Main script */
template<class Type>
Type objective_function<Type>::operator() ()
{
  Type ans=0;

  DATA_INTEGER(delay);     // Delay
  DATA_VECTOR(dt);         // Time steps
  //DATA_SCALAR(dtpredc);     // Time step for prediction
  DATA_VECTOR(dtpredcinds);
  DATA_INTEGER(dtpredcnsteps); // Number of sub time step for prediction
  //DATA_SCALAR(dtpredi);     // Time step for predictions of indices
  DATA_SCALAR(dtprediind);
  DATA_VECTOR(obsC);       // Catches
  DATA_VECTOR(ic);         // Vector such that B(ii(i)) is the state corresponding to obsC(i)
  DATA_VECTOR(nc);         // nc(i) gives the number of time intervals obsC(i) spans
  DATA_VECTOR(I);          // Index
  DATA_VECTOR(ii);         // A vector such that B(ii(i)) is the state corresponding to I(i)
  DATA_VECTOR(iq);         // A vector such that iq(i) is the index number corresponding to I_iq(i)
  DATA_VECTOR(ir);         // A vector indicating when the different rs should be used
  DATA_SCALAR(ffac);       // Management factor each year multiply the predicted F with ffac
  DATA_VECTOR(indpred);    // A vector indicating when the management factor should be applied
  DATA_SCALAR(tdfc);        // Catch Degrees of freedom of t-distribution (only used if tdf < 25)
  DATA_SCALAR(tdfi);        // Index Degrees of freedom of t-distribution (only used if tdf < 25)
  DATA_INTEGER(lamperti);  // Lamperti flag.
  DATA_INTEGER(euler);     // Euler flag.
  DATA_SCALAR(dbg);        // Debug flag, if == 1 then print stuff.
  PARAMETER(phi1);         // 
  PARAMETER(phi2);         // 
  PARAMETER(logalpha);     // sdi = alpha*sdb
  PARAMETER(logbeta);      // sdc = beta*sdf
  PARAMETER(logbkfrac);    // B0/K fraction
  PARAMETER(logF0);        // F at time 0
  PARAMETER_VECTOR(logr);  // Intrinsic growth
  PARAMETER(logK);         // Carrying capacity
  PARAMETER_VECTOR(logq);  // Catchability
  PARAMETER(logp);         // Pella-Tomlinson exponent
  PARAMETER(logsdf);       // Standard deviation for F
  PARAMETER(logsdb);       // Standard deviation for Index
  PARAMETER_VECTOR(logF);  // Random effects vector
  PARAMETER_VECTOR(logB);  // Random effects vector

  Type bkfrac = exp(logbkfrac);
  Type F0 = exp(logF0);
  int nr = logr.size();
  vector<Type> r(nr);
  for(int i=0; i<nr; i++){ r(i) = exp(logr(i)); }
  Type K = exp(logK);
  int nq = logq.size();
  vector<Type> q(nq);
  for(int i=0; i<nq; i++){ q(i) = exp(logq(i)); }
  Type p = exp(logp);
  Type sdf = exp(logsdf);
  Type sdb = exp(logsdb);
  Type sdb2 = sdb*sdb;
  Type sdi = exp(logalpha)*sdb;
  Type sdc = exp(logbeta)*sdf;
  int nobsC = obsC.size();
  int nobsCp = ic.size();
  int nobsI = I.size();
  int ns = logF.size();
  vector<Type> F = exp(logF);
  vector<Type> P(ns-1);
  vector<Type> Binf(ns);
  vector<Type> logBinf(ns);
  vector<Type> B = exp(logB);
  Type logB0 = logbkfrac + logK;
  vector<Type> Bpred(ns);
  vector<Type> rvec(ns);
  vector<Type> ffacvec(ns);
  for(int i=0; i<ns; i++){ ffacvec(i) = 1.0; }
  vector<Type> Cpred(nobsCp);
  for(int i=0; i<nobsCp; i++){ Cpred(i) = 0.0; }
  vector<Type> Cpredsub(ns);
  vector<Type> logIpred(nobsI);
  vector<Type> logCpred(nobsCp);
  vector<Type> R(nr);
  for(int i=0; i<nr; i++){ R(i) = r(i)*p/(p+1.0); }
  Type Rmean = 0.0;
  for(int i=0; i<nr; i++){ Rmean += R(i)/nr; }
  // Deterministic reference points
  Type Bmsyd = K/pow(p+1.0, 1.0/p);
  Type Fmsyd = Rmean;
  Type MSYd = Bmsyd * Fmsyd;
  Type logBmsyd = log(Bmsyd);
  Type logFmsyd = log(Fmsyd);

  // Stochastic reference points
  Type Bmsy = Bmsyd * (1.0 - (1.0 + Rmean*(p-1.0)/2.0) / (Rmean*pow(2.0-Rmean, 2.0)) * sdb2);
  Type Fmsy = Fmsyd - p*(1.0-Rmean) / pow(2.0-Rmean, 2.0) * sdb2;
  Type MSY = MSYd * (1.0 - (p+1)/2.0 / (1.0 - pow(1.0-Rmean, 2.0)) * sdb2);
  Type logBmsy = log(Bmsy);
  Type logFmsy = log(Fmsy);

  Type likval;

  if(dbg>0){
    std::cout << "--- DEBUG: script start ---" << std::endl;
    //for(int i=0; i<ns; i++) std::cout << "F(i): " << F(i) << std::endl;
    std::cout << "INPUT: logbkfrac: " << logbkfrac << std::endl;
    std::cout << "INPUT: logF0: " << logF0 << std::endl;
    for(int i=0; i<nr; i++){ std::cout << "INPUT: logr(i): " << logr(i) << " -- i: " << i << std::endl; }
    //std::cout << "INPUT: logr: " << logr << std::endl;
    std::cout << "INPUT: logK: " << logK << std::endl;
    for(int i=0; i<nq; i++){ std::cout << "INPUT: logq(i): " << logq(i) << " -- i: " << i << std::endl; }
    std::cout << "INPUT: logp: " << logp << std::endl;
    std::cout << "INPUT: logsdf: " << logsdf << std::endl;
    std::cout << "INPUT: logsdb: " << logsdb << std::endl;
    std::cout << "obsC.size(): " << obsC.size() << "  Cpred.size(): " << Cpred.size() << "  I.size(): " << I.size() << "  dt.size(): " << dt.size() << "  F.size(): " << F.size() << "  B.size(): " << B.size() << "  P.size(): " << P.size() << "  rvec.size(): " << rvec.size() << "  iq.size(): " << iq.size() << "  ic.size(): " << ic.size() << std::endl;
  }
  // Calculate rvec
  int ind;
  for(int i=0; i<ns; i++){
    ind = CppAD::Integer(ir(i)-1); // minus 1 because R starts at 1 and c++ at 0
    rvec(i) = r(ind);
    if(dbg>1){
      std::cout << "-- i: " << i << "-- ind: " << ind << " -   rvec(i): " << rvec(i) << std::endl;
    }
  }
  for(int i=1; i<indpred.size(); i++){ // don't use i=0 because this is only for plotting
    ind = CppAD::Integer(indpred(i)-1); // minus 1 because R starts at 1 and c++ at 0
    ffacvec(ind) = ffac;
    if(dbg>1){
      std::cout << "-- i: " << i << "-- ind: " << ind << " -   ffacvec(i): " << ffacvec(i) << std::endl;
    }
  }
  for(int i=0; i<ns; i++){
    if(F(i)==rvec(i)) std::cout << "Warning: F(i)-rvec(i): " << F(i)-rvec(i) << std::endl;
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
  Type sd = 1e-6; // This one is used in the hack to fix B0 and F0.


  // FISHING MORTALITY
  if(dbg>0){
    std::cout << "--- DEBUG: F loop start" << std::endl;
  }
  // Hack to set log(B(0)) equal to the fixed effect log(B0).
  likval = dnorm(logF0, logF(0), sd, 1);
  ans-=likval;
  for(int i=delay; i<ns; i++){
    Type logFpred = log(ffacvec(i)) + predictlogF(phi1, logF(i-1), phi2, logF(i-delay));
    likval = dnorm(logF(i), logFpred, sqrt(dt(i-1))*sdf, 1);
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logF(i-1): " << logF(i-1) << "  logF(i): " << logF(i) << "  ffacvec(i): " << ffacvec(i) << "  sdf: " << sdf << "  likval: " << likval << std::endl;
    }
  }

  // CALCULATE B_infinity
  for(int i=0; i<ns; i++) Binf(i) = calculateBinf(K, F(i), rvec(i), p, sdb2, lamperti); 
  logBinf = log(Binf);

  // BIOMASS PREDICTIONS
  // Hack to set log(B(0)) equal to the fixed effect log(B0).
  likval = dnorm(logB0, log(B(0)), sd, 1);
  ans-=likval;
  if(dbg>1){
    std::cout << "-- i: " << 0 << " -   logB0: " << logB0 << "  log(B(0)): " << log(B(0)) << "  sd: " << sd << "  likval: " << likval << std::endl;
  }
  for(int i=0; i<(ns-1); i++){
    // To predict B(i) use dt(i-1), which is the time interval from t_i-1 to t_i
    Bpred(i+1) = predictB(B(i), Binf(i), F(i), rvec(i), K, dt(i), p, sdb2, lamperti, euler);
    likval = dnorm(log(Bpred(i+1)), logB(i+1), sqrt(dt(i))*sdb, 1);
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logB(i+1): " << logB(i+1) << "  log(Bpred(i+1)): " << log(Bpred(i+1)) << "  sdb: " << sdb << "  likval: " << likval << std::endl;
    }
  }

  // CATCH PREDICTIONS
  for(int i=0; i<(ns-1); i++){ // ns-1 because dt is 1 shorter than state vec
    // For Cpredsub(i) use dt(i) because Cpredsub(i) is integrated over t_i to t_i+1
    Cpredsub(i) =  predictC(F(i), K, rvec(i), B(i), Binf(i), dt(i), sdb2, lamperti, euler);
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
  //if(lamperti){
  for(int i=0; i<(ns-1); i++) P(i) = B(i+1) - B(i) + Cpredsub(i);


  /*
  --- OBSERVATION EQUATIONS ---
  */

  // CATCHES
  if(dbg>0){
    std::cout << "--- DEBUG: Cpred loop start" << std::endl;
  }
  // fac and pp are used for the outlier robust Gaussian mixture.
  Type fac = 10.0;
  Type pp = 0.95;
  for(int i=0; i<nobsC; i++){
    if(tdfc < 25.0){
      //Type z = (logCpred(i)-log(obsC(i)))/sdc;
      //likval = -log(sdc) + dnorm(z, Type(0.0), Type(1.0), 1);
      //likval = -log(sdc) + ltdistr(z, Type(100.0));
      likval = pp*dnorm(logCpred(i), log(obsC(i)), sdc, 1) + (1.0-pp)*dnorm(logCpred(i), log(obsC(i)), fac*sdc, 1);
    } else {
      likval = dnorm(logCpred(i), log(obsC(i)), sdc, 1);
    }
    ans-=likval;
    // DEBUGGING
    if(dbg>1){
      std::cout << "-- i: " << i << " -   logobsC(i): " << log(obsC(i))<< "  log(Cpred(i)): " << logCpred(i) << "  sdc: " << sdc << "  likval: " << likval << std::endl;
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
      if(tdfi < 25.0){
	//Type z = (log(I(i)) - logIpred(i))/sdi;
	//likval = -log(sdi) + dnorm(z, Type(0.0), Type(1.0), 1);
	//likval = -log(sdi) + ltdistr(z, Type(100.0));
	likval = pp*dnorm(log(I(i)), logIpred(i), sdi, 1) + (1.0-pp)*dnorm(log(I(i)), logIpred(i), fac*sdi, 1);
      } else {
	likval = dnorm(log(I(i)), logIpred(i), sdi, 1);
      }
      ans-=likval;
      // DEBUGGING
      if(dbg>1){
	std::cout << "-- i: " << i << " -  ind: " << ind << " -  indq: " << indq << " -   log(I(i)): " << log(I(i)) << "  logIpred(i): " << logIpred(i) << "  sdi: " << sdi << "  likval: " << likval << std::endl;
      }
    }
  }

  // ONE-STEP-AHEAD PREDICTIONS
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
  Type logBpBmsy = logBp - logBmsy;
  Type logFp = logF(CppAD::Integer(dtprediind-1)); 
  Type logFpFmsy = logFp - logFmsy;

  vector<Type> logIp(nq);
  for(int i=0; i<nq; i++){
    logIp(i) = logq(i) + log(Bp);
  }
  
  // MSY PREDICTIONS
  /*
  Type Bpmsy;
  Type Binfpmsy;
  Type Cpmsy;
  Binfpmsy = calculateBinf(K, Fmsy, rvec(ns-1), p, sdb2, lamperti);
  Bpmsy = predictB(B(ns-1), Binfpmsy, Fmsy, rvec(ns-1), K, dtpred, p, sdb2, lamperti, euler);
  Cpmsy = predictC(Fmsy, K, rvec(ns-1), Bpmsy, Binfpmsy, dtpred, sdb2, lamperti, euler);
  Type logBpmsy = log(Bpmsy);
  */

  // Biomass and fishing mortality at last time point
  // dtpredcinds(0) is the index of B and F corresponding to the time of the last observation.
  Type logBl = logB(CppAD::Integer(dtpredcinds(0)-1)); 
  Type logBlBmsy = logBl - logBmsy;
  Type logFl = logF(CppAD::Integer(dtpredcinds(0)-1)); 
  Type logFlFmsy = logFl - logFmsy;

  // Biomass and fishing mortality over msy levels
  vector<Type> logBBmsy(ns);
  vector<Type> logFFmsy(ns);
  for(int i=0; i<ns; i++){ 
    logBBmsy(i) = logB(i) - logBmsy; 
    logFFmsy(i) = logF(i) - logFmsy; 
  }

  // ADREPORTS
  ADREPORT(r);
  ADREPORT(K);
  ADREPORT(q);
  ADREPORT(logp);
  ADREPORT(sdf);
  ADREPORT(sdc);
  ADREPORT(sdb);
  ADREPORT(sdi);
  ADREPORT(MSY);
  ADREPORT(MSYd);
  // B
  ADREPORT(Bmsy);
  ADREPORT(logBmsy);
  ADREPORT(Bmsyd);
  ADREPORT(logBmsyd);
  ADREPORT(logBp);
  //ADREPORT(logBpmsy);
  ADREPORT(logBpBmsy);
  ADREPORT(logBinf);
  ADREPORT(logB0);
  ADREPORT(logBl);
  ADREPORT(logBlBmsy);
  ADREPORT(logBBmsy);
  // F
  ADREPORT(Fmsy);
  ADREPORT(logFmsy);
  ADREPORT(Fmsyd);
  ADREPORT(logFmsyd);
  ADREPORT(logFp);
  ADREPORT(logFpFmsy);
  ADREPORT(logFl);
  ADREPORT(logFlFmsy);
  ADREPORT(logFFmsy);
  // C
  //ADREPORT(Cpmsy);
  ADREPORT(Cpredsub);
  ADREPORT(logCpred);
  ADREPORT(logCp);
  // Other
  ADREPORT(logIpred);
  ADREPORT(P);
  // REPORTS (these don't require sdreport to be output)
  REPORT(Cp);
  REPORT(logIp);

  return ans;
}

