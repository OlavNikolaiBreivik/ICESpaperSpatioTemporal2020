#' Run model
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb
#' @param data Data
#' @param par Parameters
#' @param conf Configurations
#' @param doDetailedRep Return reported indices etc if true.
#' @param low Bounds
#' @param up Bounds
#' @useDynLib spatioTemporalIndices
#' @return A fitted stim object
#' @details
#' @export
fitModel<-function(data,par,conf,map = NULL, doDetailedRep = 2, low = list(),up = list()){

  tryCatch(
    {
      #Will fail if not using MKL (efficient linear algebra library)
      setMKLthreads(1) #not profiting much by using more cores
    },
    error=function(cond) {
      message("MKL library not used, couputation will be time consuming")
      return(NA)
    }
  )

  #Set map
  if(is.null(map)){
    map = setMap(par,conf)
  }

  #Configurations used on C-side
  data$spatial=conf$spatial
  data$spatioTemporal=conf$spatioTemporal
  data$splineDepth=conf$splineDepth[2]
  data$useNugget=conf$nugget
  data$pcPriorsRange = conf$pcPriorRange
  data$pcPriorsSD = conf$pcPriorsd
  data$usePCpriors = conf$usePcPriors
  data$zeroInflated = conf$zeroInflated
  data$doDetailedRep = doDetailedRep



  obj <- MakeADFun(data, par, random=c("xS","xST","betaDepth", "nugget","nuggetIndex"),profile = c("beta0","betaSun"), DLL="spatioTemporalIndices",map = map)

  lower<-rep(-Inf,length(obj$par))
  upper<-rep(Inf,length(obj$par))
  for(nn in names(low)) lower[names(obj$par)==nn]=low[[nn]]
  for(nn in names(up)) upper[names(obj$par)==nn]=up[[nn]]

  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(trace = 1,iter.max = 1000, eval.max = 1000),
                lower = lower, upper = upper)


  if(doDetailedRep){
    rep <- sdreport(obj)
  }else{
    rep <- sdreport(obj,ignore.parm.uncertainty  = TRUE)
  }

  toReturn = list(obj = obj,opt = opt,rep = rep,conf = conf,data = data,map = map,par = par)
  class(toReturn) = "stim"

  return(toReturn)
}



#' jit
#' @param run The result of running fitModel
#' @param njit Number of jitter runs
#' @param ncores Number of cores to use
#' @param sd Standard deviation of noise to start values
#' @details
#' @return
#' @export
#'
jit<-function(run,njit,ncores = 1,sd = 0.2){


  #Construct jitter
  parOriginal = setPar(run$data,run$conf)
  pOriginal <- unlist(parOriginal)
  par <- lapply(1:njit, function(i)relist(pOriginal+rnorm(length(pOriginal),sd=sd), run$par))

  #Set to spatio-temporal parameters to arrays as needed in the current implementation
  for(i in 1:length(par)){
    par[[i]]$xS = array(par[[i]]$xS, dim = dim(parOriginal$xS))
    par[[i]]$xST = array(par[[i]]$xST, dim = dim(parOriginal$xST))
  }

  #Set those who are mapped as NA to its initial value
  par = lapply(par, function(p){
    for(i in 1:length(run$map)){
      j = which(names(p) == names(run$map)[i])
      p[[j]][which(is.na(run$map[[i]]))] = parOriginal[[j]][which(is.na(run$map[[i]]))]
    }
    p
  })

  if(ncores ==1){
    runs=lapply(par,function(f) fitModel(data = run$data,par = f,conf = run$conf, map = run$map))
  }else{
    cl <-makeCluster(min(njit,ncores),outfile = "")
    clusterExport(cl,varlist=c("run","par"),envir=environment())
    clusterEvalQ(cl, library("spatioTemporalIndices"))
    runs=parLapply(cl,par,function(f) fitModel(data = run$data,par = f,conf = run$conf, map = run$map))
    stopCluster(cl)
  }

  attributes(runs)$runOriginal = run

  p <- lapply(runs, function(f){
    unlist(as.list(f$rep,"Est"))
  })

  #Max diff of indices and log-likelihood
  rl = as.list(run$rep, "Est", report = TRUE)

  maxDiffIndices = lapply(runs,function(f){
    rlJ = as.list(f$rep, "Est", report = TRUE)
    max(abs(rlJ$muReportFull- rl$muReportFull))
  })

  maxDiffLogLik = lapply(runs,function(f){
    max(abs(f$opt$objective- run$opt$objective))
  })

  maxDiffSun = lapply(runs,function(f){
    rlJ = as.list(f$rep, "Est", report = TRUE)
    max(abs(rlJ$fourierReportLow- rl$fourierReportLow))
  })

  maxDiffDepth= lapply(runs,function(f){
    rlJ = as.list(f$rep, "Est", report = TRUE)
    max(abs(rlJ$depthReport1- rl$depthReport1))
  })

  pl = as.list(run$rep,"Est")
  pl$xS = pl$xS* exp(pl$log_sigma[1])^2
  pl$xST = pl$xST* exp(pl$log_sigma[2])^2
  pl$nugget = pl$xST* exp(pl$log_sigma[3])^2


  diffPar = lapply(runs, function(f) {
    plJit = as.list(f$rep, "Est")
    plJit$xS = plJit$xS* exp(plJit$log_sigma[1])^2
    plJit$xST = plJit$xST* exp(plJit$log_sigma[2])^2
    plJit$nugget = plJit$xST* exp(plJit$log_sigma[3])^2
    mpd = mapply(function(g,h){
      max(abs(g-h))
    } , plJit, pl)
  })

  maxMat = t(matrix(unlist(diffPar), nrow = length(diffPar[[1]])))
  maxVec = apply(maxMat,2, max)

  maxVecAll = c(maxVec, max(unlist(maxDiffIndices)),
                    max(unlist(maxDiffLogLik)),
                    max(unlist(maxDiffSun)),
                    max(unlist(maxDiffDepth)))
  names(maxVecAll) = c(names(pl), "Index", "logLik", "Sun","Depth")


  return(list(runs = runs, maxVecAll = maxVecAll,maxMat = maxMat))
}


#' runwithout Fit stim and leave out observation indices as provided in fold
#' @param data Data
#' @param par Parameters
#' @param conf Configurations
#' @param map map-argument
#' @param fold Vector of data indices to remove when fitting stim
#' @details
#' @export
runwithout = function(data,par,conf,map,fold){
  predMatrix = data$fishObsMatrix*0
  predMatrix[fold,] = 1

  data$predMatrix = predMatrix
  run = fitModel(data = data,par = par,conf = conf,map=map)
  return(run)
}



