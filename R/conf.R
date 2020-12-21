
##' defConf
##'
##' Configurations used. The object returned is
##' used in \code{\link{setupData}} and in  \code{\link{setPar}}.
##'
##' @param years Years of data included
##' @param lengthGroups Length groups included
##' @param spatial If 0-no spatial effect, 1-Include spatial effect
##' @param spatioTemporal If 0-no spatial-temporal effect, 1-Include spatial-temporal effect
##' @param nugget If 0-no nugget effect, 1-Include nugget effect
##' @param beta0 0-no intercept (probably not used), 1- intercept included
##' @param splineDepth Vector with depth configurations, first element is the k-variable in spline (6 means 6 basis function). Second element: 0- no depth effect, 1-depth effectk, 2-length dependent depth effect.
##' @param sunAlt Vector with sun altitude configurations, first element is the numberof basis functions. Second element: 0-sun effectk, 1-length dependent sun effect.
##' @param duration Duration effect, 0- used as offset, 1- used as covariate.
##' @param maxLength Maximum length, used when defining the pluss group
##' @param minLength Minimum length, used when defining the minus group (probably never used)
##' @param reduceLength The resolution in length dimension used in latent effect
##' @param cutoff Cutoff used in mesh, first element is for spatial effect, second for spatio-temporal effect, third is 0 if same mesh is used for both spatial and spatio-temporal effect
##' @param cbound cbound used in mesh, first element is for spatial effect, second for spatio-temporal effect
##' @param pcPriorRange pc-priors for spatial and spatio-temporal effect
##' @param pcPriorsd pc-priors for spatial and spatio-temporal effect
##' @param usePcPriors If 0-No pc-priors used, 1- use pc-priors (probably never used)
##' @param zeroInflated If 0-No zero inflation, 1- use zero inflation (NB! not fully implemented)
##' @details
##' @export
defConf <- function(years, lengthGroups,spatial = 1,spatioTemporal = 1,nugget = 1,beta0=1,splineDepth=c(6,1),sunAlt=c(1,0),duration=0,
                    maxLength = 100,minLength = 5,reduceLength = 3,
                    cutoff = c(100,100,0), cbound = c(400,100),
                    pcPriorRange = c(100,0.5),pcPriorsd = c(0.5,0.5), usePcPriors = 0, zeroInflated = 0,
                    obsModel = 2,
                    mapRhoL = c(0,1,2), simulateProcedure = 1){
  conf= list()
  conf$lengthGroups = lengthGroups
  conf$years = years
  conf$minLength=minLength
  conf$maxLength=maxLength
  conf$dLength = 5

  #Covariates
  conf$sunAlt = sunAlt
  conf$duration =duration
  conf$splineDepth=splineDepth
  conf$beta0=beta0
  conf$spatial = spatial
  conf$spatioTemporal = spatioTemporal
  conf$cutoff = cutoff
  conf$cbound = cbound

  lTmp = rep(1:100,each = reduceLength)
  conf$lengthGroupsReduced = lTmp[1:length(conf$lengthGroups)]
  conf$reduceLength = reduceLength

  conf$pcPriorRange = pcPriorRange
  conf$pcPriorsd= pcPriorsd
  conf$usePcPriors = usePcPriors

  conf$zeroInflated = zeroInflated

  conf$obsModel = obsModel
  conf$nugget = nugget
  conf$mapRhoL = mapRhoL
  conf$simulateProcedure = simulateProcedure
  return(conf)
}

##' defConfPred
##'
##'Configurations used for prediction.
##'
##'
##' @param yearToPredict Prediciton year
##' @param strata Which statas to predict indices in (probably always all stratas and extract more detailed predictions on a later stage)
##' @param conf Configurations used when fitting the model
##' @param n provide number of integration points and number of bootstrap samples in internal validation routine
##' @details
##' @export
defConfPred <- function(strata = 1:26,nIntPoints= 4000){
  confPred = list()
  confPred$Strata=strata
  confPred$nIntPoints = nIntPoints
  return(confPred)
}

##' setMap
##'
##'\code{setMap} defines the map-argument used in MakeADFun.
##'
##' @param par Parameters included
##' @param conf Configurations
##' @details
##' @export
setMap <- function(par, conf){
  map= list()
  map$log_kappa = as.factor(c(0,1))
  map$tan_rho_l = conf$mapRhoL #Default use same parameter for length correlation

  if(conf$nugget[1] ==0){
    map$log_sigma = as.factor(c(0,1,NA))
    map$nugget = as.factor(rep(NA,length(par$nugget)))
    map$nuggetIndex = as.factor(rep(NA,length(par$nuggetIndex)))
    map$tan_rho_l[3] = NA
  }else{
    map$log_sigma = as.factor(c(0,1,2))
  }


  if(conf$sunAlt[1]==0){
    map$betaSun = as.factor(rep(NA,length(par$betaSun))) #Not use time in day
  }else if(conf$sunAlt[2]==0){
    tmp = 0:(conf$sunAlt[1]*2-1)
    map$betaSun = as.factor(c(tmp,tmp)) #Not use length dependent time in day
  }
  if(conf$duration==0){
    map$betahaul = as.factor(rep(NA,length = length(par$betahaul)))
  }
  if(length(conf$years)==1){
    map$tan_rho_t = as.factor(NA)
    if(conf$spatial==1 & conf$spatioTemporal==1){
      warning("Using both spatial and spatio-temporal with only one year of data will overparametrize the model, turn of spatio-temporal contribution")
      conf$spatioTemporal=0;
    }
  }
  if(conf$splineDepth[2]==0){
    map$betaDepth=as.factor(rep(NA,length(par$betaDepth)))
    map$log_lambda=as.factor(c(NA,NA))
  }else if(conf$splineDepth[2]==1){
    tmp = 0:(length(par$betaDepth)/2-1)
    map$log_lambda=as.factor(c(0,0))
    map$betaDepth = as.factor(c(tmp,tmp)) #Not use length dependent depth effect
  }else if(conf$splineDepth[2]==1){
    map$log_lambda=as.factor(c(0,0))#Use same lambda in both depth splines
  }
  if(conf$beta0==0){
    map$beta0=par$beta0
    for(i in 1:dim(map$beta0)[2]){
      map$beta0[,i]=i
    }
    map$beta0=as.factor(map$beta0)
  }

  if(conf$spatial ==0){
    map$xS = as.factor(rep(NA,length(par$xS)))
    if(conf$nugget[1] ==0){
      map$log_sigma= as.factor(c(NA,0,NA))
    }else{
      map$log_sigma= as.factor(c(NA,0,1))
    }
    map$log_kappa= as.factor(c(NA,0))
    map$tan_rho_l[1]= NA
    map$tan_rho_l[2:3]  = map$tan_rho_l[2:3]-1
  }
  if(conf$spatioTemporal ==0){
    map$xST = as.factor(rep(NA,length(par$xST)))
    map$tan_rho_t = as.factor(NA)
    if(conf$spatial ==0){
      map$log_kappa= as.factor(c(NA,NA))
      if(conf$nugget[1] ==0){
        map$log_sigma= as.factor(c(NA,NA,NA))
        map$tan_rho_l[3]= NA
      }else{
        map$log_sigma= as.factor(c(NA,NA,1))
        map$tan_rho_l[3]= 0
      }
      map$tan_rho_l[1:2]= c(NA,NA)

    }else{
      map$log_kappa= as.factor(c(0,NA))
      if(conf$nugget[1] ==0){
        map$log_sigma= as.factor(c(0,NA,NA))
        map$tan_rho_l[3] =   NA
      }else{
        map$log_sigma= as.factor(c(0,NA,1))
        map$tan_rho_l[3] =   map$tan_rho_l[3]-1
      }
      map$tan_rho_l[2] =  NA
    }
  }

  if(conf$obsModel==2){
    map$logSize = as.factor(NA)
  }



  if(conf$zeroInflated==0){
    map$delta_z = as.factor(rep(NA,length(par$delta_z)))
  }else{
    warning("Map-functionality not yet implemented for selected zero-inflation procedure.")
  }

  map$tan_rho_l = as.factor(map$tan_rho_l)
  return(map)
}

