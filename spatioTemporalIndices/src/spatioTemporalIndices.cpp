#define TMB_LIB_INIT R_init_spatioTemporalIndices

#include <TMB.hpp>
#include <cmath>
using namespace tmbutils;


/* List of sparse matrices */
template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      if(!isValidSparseMatrix(sm))
            error("Not a sparse matrix");
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density;

  //Load data and parameters
  DATA_MATRIX(fishObsMatrix); //Observations in a matrix
  DATA_VECTOR(dist); // Distance of each haul
  DATA_STRUCT(A_ListS, LOSM_t);
  DATA_STRUCT(A_ListST, LOSM_t);
  DATA_STRUCT(spdeMatricesS,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_STRUCT(spdeMatricesST,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_VECTOR(nStationsEachYear)
  DATA_INTEGER(numberOfLengthGroups);
  DATA_SPARSE_MATRIX(S_depth);
  DATA_INTEGER(Sdim);
  DATA_VECTOR(lengthGroupsReduced);//used when reducing length dimension
  DATA_MATRIX(predMatrix); //used to leave out data from likelihood
  DATA_MATRIX(X_sunAlt);
  DATA_MATRIX(X_sunAltReport);
  DATA_MATRIX(X_sunAltIntegrate);
  DATA_MATRIX(X_depth);
  DATA_MATRIX(X_depthReport);
  DATA_MATRIX(X_depth_int);
  DATA_VECTOR(weigthLength); //wieght used when reducing length dimension
  DATA_VECTOR(areas); //Areas of stratas
  DATA_VECTOR(selectedStratas); //Stratas to calculate the joint index from
  DATA_INTEGER(zeroInflated); //if !=0, zero inflated model is used, work in progress
  DATA_INTEGER(nBasisSunAlt);
  DATA_INTEGER(obsModel);
  DATA_INTEGER(doDetailedRep);
  DATA_INTEGER(doValidation);


  DATA_VECTOR(obsVector);//Observations in a vector. TODO: remove the observation matrix.
  DATA_IVECTOR(idxStart);//Index helper when refering to the observation vector
  DATA_VECTOR_INDICATOR(keep,obsVector);

  DATA_VECTOR(xInt); //Used to calcualte COF
  DATA_VECTOR(yInt); //Used to calcualte COF

  DATA_SPARSE_MATRIX(ApredS);
  DATA_SPARSE_MATRIX(ApredST);
  DATA_INTEGER(nStrata);
  DATA_MATRIX(predLoc);

  DATA_INTEGER(spatial); //Configuration spatial, include if 1
  DATA_INTEGER(spatioTemporal); //Configuration spatio-temporal, include if 1
  DATA_INTEGER(splineDepth); //Configuration spatio-temporal, include if 1
  DATA_INTEGER(useNugget); //Configuration useNugget, include if 1

  DATA_INTEGER(simulateProcedure);

  //Pc-priors for latent effect
  DATA_VECTOR(pcPriorsRange);
  DATA_VECTOR(pcPriorsSD);
  DATA_INTEGER(usePCpriors);

  PARAMETER_MATRIX(beta0);
  PARAMETER_VECTOR(betahaul);
  PARAMETER_VECTOR(betaSun);
  PARAMETER_VECTOR(betaDepth);
  PARAMETER_VECTOR(log_lambda);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(log_kappa);
  PARAMETER(logSize);
  PARAMETER_VECTOR(tan_rho_t);
  PARAMETER_VECTOR(tan_rho_l);
  PARAMETER_VECTOR(delta_z);

  PARAMETER_ARRAY(xS);
  PARAMETER_ARRAY(xST);
  PARAMETER_ARRAY(nugget);
  PARAMETER_ARRAY(nuggetIndex);


  //Transform parameters
  vector<Type> sigma = exp(log_sigma);
  vector<Type> kappa = exp(log_kappa);
  Type size= exp(logSize);

  vector<Type> rho_t=tan_rho_t;
  rho_t(0)=Type(2)/(Type(1) + exp(-Type(2) * tan_rho_t(0))) - Type(1);

  vector<Type> rho_l=tan_rho_l;
  rho_l(0) = Type(2)/(Type(1) + exp(-Type(2) * tan_rho_l(0))) - Type(1);
  rho_l(1) = Type(2)/(Type(1) + exp(-Type(2) * tan_rho_l(1))) - Type(1);
  rho_l(2) = Type(2)/(Type(1) + exp(-Type(2) * tan_rho_l(2))) - Type(1);

  vector<Type> lambda = exp(log_lambda);

  for(int i =1; i<delta_z.size(); ++i){
    delta_z(i) = exp(delta_z(i));
  }

  //Construct Q in spatial dimension
  SparseMatrix<Type> Q_s = Q_spde(spdeMatricesS,kappa(0));
  SparseMatrix<Type> Q_st = Q_spde(spdeMatricesST,kappa(1));

  //Calculates nll-------------------------------
  Type nll = 0.0;
 // parallel_accumulator<Type> nll(this);

  //Latent effects
  Type d = 2; //Part of spatial pc-prior
  Type rhoP;
  Type R = -log(pcPriorsRange(1))*pow(pcPriorsRange(0),d/2);
  Type S = -log(pcPriorsSD(1))/pcPriorsSD(0);
  if(spatial==1){
    if(usePCpriors==1){
      rhoP = sqrt(8)/kappa(0);
      nll -= log( d/2 * R *S * pow(rhoP,(-1-d/2))* exp(-R* pow(rhoP,(-d/2)) -S* sigma(0)) ); //pc-prior contribution
    }
    nll += SEPARABLE(AR1(rho_l(0)),GMRF(Q_s))(xS);
    if(simulateProcedure==1){
      SIMULATE{
        SEPARABLE(AR1(rho_l(0)),GMRF(Q_s)).simulate(xS);
      }
    }
  }
  if(spatioTemporal==1){
    if(usePCpriors==1){
      rhoP = sqrt(8)/kappa(1);
      nll -= log( d/2 * R *S * pow(rhoP,(-1-d/2))* exp(-R* pow(rhoP,(-d/2)) -S* sigma(1)) ); //pc-prior contribution
    }
    nll += SEPARABLE(AR1(rho_l(1)),SEPARABLE(AR1(rho_t(0)),GMRF(Q_st)))(xST);
    if(simulateProcedure==1){
      SIMULATE{
        SEPARABLE(AR1(rho_l(1)),SEPARABLE(AR1(rho_t(0)),GMRF(Q_st))).simulate(xST);
      }
    }
  }

  if(useNugget==1){
    int nHaul = CppAD::Integer(nStationsEachYear.sum());

    SparseMatrix<Type> Q_nuggetIID(nHaul,nHaul);
    for(int i = 0; i< nHaul; ++i){
      Q_nuggetIID.coeffRef(i,i)=1;
    }

    nll += SEPARABLE(GMRF(Q_nuggetIID),AR1(rho_l(2)))(nugget);
    if(simulateProcedure==1){
      SIMULATE{
        SEPARABLE(GMRF(Q_nuggetIID),AR1(rho_l(2))).simulate(nugget);
      }
    }
    SparseMatrix<Type> Q_nuggetIIDI(xInt.size(),xInt.size());
    for(int i = 0; i< xInt.size(); ++i){
      Q_nuggetIIDI.coeffRef(i,i)=1;
    }

    nll += SEPARABLE(GMRF(Q_nuggetIIDI),AR1(rho_l(2)))(nuggetIndex);

  }

  //p-spline
  int nSplineDepth = betaDepth.size()/2;
  vector<Type> parDepth1(nSplineDepth);
  vector<Type> parDepth2(nSplineDepth);
  for(int i =0; i<nSplineDepth; ++i){
    parDepth1(i) = betaDepth(i);
    parDepth2(i) = betaDepth(i + nSplineDepth);
  }
  if(splineDepth ==1){
    nll -= Type(0.5)*Sdim*log_lambda(0) - 0.5*lambda(0)*GMRF(S_depth).Quadform(parDepth1);
  }else if(splineDepth ==2){ //Two length dependent splines
    nll -= Type(0.5)*Sdim*log_lambda(0) - 0.5*lambda(0)*GMRF(S_depth).Quadform(parDepth1);
    nll -= Type(0.5)*Sdim*log_lambda(1) - 0.5*lambda(1)*GMRF(S_depth).Quadform(parDepth2);
  }

  //fourier
  vector<Type> betaSunLow(nBasisSunAlt*2);
  vector<Type> betaSunHigh(nBasisSunAlt*2);
  for(int i =0; i<(nBasisSunAlt*2); ++i){
    betaSunLow(i) = betaSun(i);
    betaSunHigh(i) = betaSun(i + nBasisSunAlt*2);
  }


  Type log_var_minus_mu;
  Type covariatesConvexW; //Coefficient in convex combination of length dependent effects
  matrix<Type> mu(fishObsMatrix.rows(),fishObsMatrix.cols());
  vector<Type> deltaS, deltaS2, deltaST, deltaST2; //Smoothing in length dimension
  matrix<Type> deltaMatrixS(numberOfLengthGroups,999), deltaMatrixST(numberOfLengthGroups,999);//Smoothed effect in length dimension
  vector<Type> timeInDayLow = X_sunAlt*betaSunLow;
  vector<Type> timeInDayHigh = X_sunAlt*betaSunHigh;
  vector<Type> depthEffect1=X_depth*parDepth1;
  vector<Type> depthEffect2=X_depth*parDepth2;

  vector<Type> validation(numberOfLengthGroups);
  validation.setZero();
  int counter=0;
  for(int y=0;y<nStationsEachYear.size();y++){
    SparseMatrix<Type> As = A_ListS(y);
    SparseMatrix<Type> Ast = A_ListST(y);
    //Latent effects contribution
    if(CppAD::Integer(lengthGroupsReduced(0))==CppAD::Integer(lengthGroupsReduced(1))){ //Length dimension is reduced
      for(int l=0; l <numberOfLengthGroups;++l){
        deltaS = As * xS.col(CppAD::Integer(lengthGroupsReduced(l))).matrix();
        deltaS2 = As * xS.col(CppAD::Integer(lengthGroupsReduced(l))+1).matrix();
        deltaST = Ast * xST.col(CppAD::Integer(lengthGroupsReduced(l))).col(y).matrix();
        deltaST2 = Ast * xST.col(CppAD::Integer(lengthGroupsReduced(l))+1).col(y).matrix();

        for(int s=0; s<nStationsEachYear(y);++s){
          deltaMatrixS(l,s) = weigthLength(l)*deltaS(s) + (1-weigthLength(l))*deltaS2(s);
          deltaMatrixST(l,s) = weigthLength(l)*deltaST(s) + (1-weigthLength(l))*deltaST2(s);
        }
      }
    }else{
      for(int l=0; l <numberOfLengthGroups;++l){
        deltaS = As * xS.col(CppAD::Integer(lengthGroupsReduced(l))).matrix();
        deltaST = Ast * xST.col(CppAD::Integer(lengthGroupsReduced(l))).col(y).matrix();
        for(int s=0; s<nStationsEachYear(y);++s){
          deltaMatrixS(l,s) = deltaS(s);
          deltaMatrixST(l,s) = deltaST(s);
        }
      }
    }

    Type muZero;
    Type pZero;
    Type pPos;
    for(int s=0; s<nStationsEachYear(y);++s){
      for(int l=0; l <numberOfLengthGroups;++l){
	         covariatesConvexW = (numberOfLengthGroups-l-1)/(numberOfLengthGroups-1);
           mu(counter,l)= exp( beta0.row(y)(l)+
             covariatesConvexW*timeInDayLow(counter) + (1-covariatesConvexW)*timeInDayHigh(counter) +
             deltaMatrixS.row(l)(s)*sigma(0)*sigma(0)+
             deltaMatrixST.row(l)(s)*sigma(1)*sigma(1)+
             betahaul(l)*log(dist(counter))+
             covariatesConvexW*depthEffect1(counter) + (1-covariatesConvexW)*depthEffect2(counter)+
             nugget.col(counter)(l)*sigma(2)*sigma(2));
        log_var_minus_mu=log(mu(counter,l)*mu(counter,l)*size);
        if(predMatrix(counter,l)==0){
	        if(zeroInflated ==1){
	          muZero = exp(delta_z(0) +
	            delta_z(1)* beta0.row(y)(l) +
	            delta_z(2)* (covariatesConvexW*timeInDayLow(counter) + (1-covariatesConvexW)*timeInDayHigh(counter)) +
	            delta_z(3)* deltaMatrixS.row(l)(s)*sigma(0)*sigma(0)+
	            delta_z(4)* deltaMatrixST.row(l)(s)*sigma(1)*sigma(1)+
	            delta_z(5)* (covariatesConvexW*depthEffect1(counter) + (1-covariatesConvexW)*depthEffect2(counter))+
	            delta_z(7)* nugget.col(counter)(l)*sigma(2)*sigma(2)+
	            delta_z(8)* betahaul(l)*log(dist(counter)));
	          pZero = dpois(Type(0), muZero,true);
	          if(obsModel==1){
	            pPos = dnbinom_robust(obsVector(idxStart(counter) +l), log(mu(counter,l)),log_var_minus_mu,true)  + logspace_sub(Type(0),pZero);
	          }else if (obsModel==2){
	            pPos = dpois(obsVector(idxStart(counter) +l), mu(counter,l),true)  + logspace_sub(Type(0),pZero);
	          }else{
	            //Not implemented
	            exit(0);
	          }

	          if(fishObsMatrix(counter,l)==0){
	            nll -=keep(idxStart(counter) +l)*logspace_add(pZero,pPos);
	          }else{
	            nll -=keep(idxStart(counter) +l)*pPos;
	          }
	        }else{
	          if(obsModel==1){
	            nll -= keep(idxStart(counter) +l)*dnbinom_robust(obsVector(idxStart(counter) +l), log(mu(counter,l)),log_var_minus_mu,true);
	            SIMULATE{
	              fishObsMatrix(counter,l) = rnbinom2(mu(counter,l),mu(counter,l) + mu(counter,l)*mu(counter,l)/size);
	            }
	          }else if(obsModel==2){
              nll -= keep(idxStart(counter) +l)*dpois(obsVector(idxStart(counter) +l), mu(counter,l),true);
	            SIMULATE{
	              fishObsMatrix(counter,l) = rpois(mu(counter,l));
	            }
	          }
	        }
        }else{
          validation(l) = validation(l) + mu(counter,l);
        }
    	}
      counter++;
	  }
  }
  //---------------------------------------------

  if(simulateProcedure==1){
    SIMULATE{
      REPORT(fishObsMatrix);
      REPORT(nugget);
      REPORT(xS);
      REPORT(xST);
    }
  }else if(simulateProcedure==2){
    REPORT(fishObsMatrix);
  }

  if(doValidation==1){
    ADREPORT(validation);
  }

  //Calculate indecies and COG
  int nYears = nStationsEachYear.size();
  array<Type>  muReport(nYears,nStrata,fishObsMatrix.cols()); //Log-index in each strata
  matrix<Type> muReportFull(nYears,fishObsMatrix.cols());     //Log-index in survey area
  matrix<Type> muReportSelected(nYears,fishObsMatrix.cols()); //Log-index in selected stratas
  matrix<Type> COGx(nYears,fishObsMatrix.cols());
  matrix<Type> COGy(nYears,fishObsMatrix.cols());
  matrix<Type> muReportSelectedExp(nYears,fishObsMatrix.cols());//Index in survey area
  matrix<Type> muReportFullExp(nYears,fishObsMatrix.cols()); //Index in selected stratas

  vector<Type> deltaPredST;
  vector<Type> deltaPredS;

  COGx.setZero();COGy.setZero();muReportFull.setZero();muReportSelected.setZero();muReport.setZero();

  vector<Type> depthEffectInt1=X_depth_int*parDepth1;
  vector<Type> depthEffectInt2=X_depth_int*parDepth2;
  vector<Type> timeInDayEffectIntLow = X_sunAltIntegrate*betaSunLow;
  vector<Type> timeInDayEffectIntHigh = X_sunAltIntegrate*betaSunHigh;

  for(int y=0; y<nYears; ++y){
    for(int l =0; l<numberOfLengthGroups; ++l){
      covariatesConvexW = (numberOfLengthGroups-l-1)/(numberOfLengthGroups-1);
      if(CppAD::Integer(lengthGroupsReduced(0))==CppAD::Integer(lengthGroupsReduced(1))){
        deltaPredS = weigthLength(l)*ApredS * xS.col(CppAD::Integer(lengthGroupsReduced(l))).matrix()+
          (1-weigthLength(l))*ApredS * xS.col(CppAD::Integer(lengthGroupsReduced(l))+1).matrix();
        deltaPredST = weigthLength(l)*ApredST * xST.col(CppAD::Integer(lengthGroupsReduced(l))).col(y).matrix()+
          (1-weigthLength(l))*ApredST * xST.col(CppAD::Integer(lengthGroupsReduced(l))+1).col(y).matrix();
      }else{
        deltaPredS = ApredS * xS.col(CppAD::Integer(lengthGroupsReduced(l))).matrix();
        deltaPredST = ApredST * xST.col(CppAD::Integer(lengthGroupsReduced(l))).col(y).matrix();
      }
      Type COGtotal = 0;
      Type muThis = 0;
      for(int strata =1; strata<(nStrata+1); ++strata){
        vector<Type> predObsInStrata = predLoc.row(strata-1);
        Type nIntegrate = 0;
        muReport(y,strata-1,l) = 0;
        for(int i =0; i<predObsInStrata.size(); ++i){
          if(predObsInStrata(i)>=0){
            muThis =  exp(beta0.row(y)(l) +
              covariatesConvexW*timeInDayEffectIntLow(0) + (1-covariatesConvexW)*timeInDayEffectIntHigh(0)+
              covariatesConvexW*depthEffectInt1(CppAD::Integer(predObsInStrata(i))) + (1-covariatesConvexW)*depthEffectInt2(CppAD::Integer(predObsInStrata(i)))+
              deltaPredS(CppAD::Integer(predObsInStrata(i)))*sigma(0)*sigma(0)+
              deltaPredST(CppAD::Integer(predObsInStrata(i)))*sigma(1)*sigma(1)+
              nuggetIndex.col(CppAD::Integer(predObsInStrata(i)))(l)*sigma(2)*sigma(2) //Should replace with 0.5sigma^2
              );

            muReport(y,strata-1,l)  =muReport(y,strata-1,l) +muThis;

            COGx(y,l) +=  muThis *xInt(CppAD::Integer(predObsInStrata(i)));
            COGy(y,l) +=  muThis *yInt(CppAD::Integer(predObsInStrata(i)));
            COGtotal += muThis;

            nIntegrate += 1;
          }
        }
        muReportFull(y,l) += muReport(y,strata-1,l) *areas(strata-1)/nIntegrate;

        bool  insideSelected = false;
        for(int ss = 0; ss<selectedStratas.size(); ++ss){
          if(CppAD::Integer(selectedStratas(ss)) == strata){
            insideSelected= true;
          }
        }
        if(insideSelected){
          muReportSelected(y,l) += muReport(y,strata-1,l) *areas(strata-1)/nIntegrate;
        }

        muReport(y,strata-1,l)  = log(muReport(y,strata-1,l) *areas(strata-1)/nIntegrate);
      }
      muReportFullExp(y,l) = muReportFull(y,l);
      muReportSelectedExp(y,l) = muReportSelected(y,l);
      muReportFull(y,l) = log(muReportFull(y,l));
      muReportSelected(y,l) = log(muReportSelected(y,l));

      COGx(y,l) = COGx(y,l)/COGtotal;
      COGy(y,l) = COGy(y,l)/COGtotal;
    }
  }



  if(doDetailedRep==1){
    ADREPORT(muReportFullExp); //Indices which are bias corrected
    ADREPORT(muReportSelectedExp); //Indices which are bias corrected

    //ADREPORT(muReport); //Report all stratas, remove for computational reasons
    ADREPORT(muReportFull);
    ADREPORT(muReportSelected);

    ADREPORT(COGx);
    ADREPORT(COGy);

  }else if(doDetailedRep==2){
    ADREPORT(muReportFullExp); //Indices which are bias corrected
    ADREPORT(muReportSelectedExp); //Indices which are bias corrected

    //ADREPORT(muReport); //Report all stratas, remove for computational reasons
    ADREPORT(muReportFull);
    ADREPORT(muReportSelected);

    matrix<Type> COGne = COGx + COGy;
    ADREPORT(COGx);
    ADREPORT(COGy);

    vector<Type> fourierReportLow = X_sunAltReport*betaSunLow;
    vector<Type> fourierReportHigh = X_sunAltReport*betaSunHigh;
    ADREPORT(fourierReportLow);
    ADREPORT(fourierReportHigh);

    vector<Type> depthReport1 =X_depthReport*parDepth1;
    vector<Type> depthReport2 =X_depthReport*parDepth2;
    ADREPORT(depthReport1);
    ADREPORT(depthReport2);
  }
  return nll;
}
