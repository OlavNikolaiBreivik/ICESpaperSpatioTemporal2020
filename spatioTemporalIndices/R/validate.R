#' validateModel
#' @param
#' @return
#' @export
#' @examples
#' @return Validation suggested by reviewer2
validateModel<-function(run,nCores = 1){
  dataDir <- system.file("extdata", package = "spatioTemporalIndices")
  preStratas=read.delim(paste0(dataDir,"/vintertokt_barentshavny.txt"),head=F)

  #Reviewer 2's suggestion.
  specifiedStratas=c(7,8,9,10,13,14,15,16,17,20)

  kFold = list()
  counter = 1
  for(y in run$conf$years){
      indexY = which(attributes(run$data)$year==y)
      locObsLatLonY =attributes(run$data)$locObsLatLon[indexY,]
      kFold[[counter]] = pointsInStratum(strataNumbers = specifiedStratas,stratas = preStratas,locations=locObsLatLonY) + min(indexY)-1
      counter = counter + 1
  }

  if(nCores==1){
    runs <- lapply(kFold, function(f)runwithout(run$data,run$par,run$conf,run$map, fold=f,doDetailedRep = 0, doValidation = 1))
  }else{
    print(paste0("use: ",nCores, " cores"))
    data = run$data
    pl = run$pl
    conf = run$conf
    map = run$map
    cl <-makeCluster(nCores,outfile = "")
    clusterExport(cl,varlist=c("data","pl","conf","map", "kFold"),envir=environment())
    clusterEvalQ(cl, library("spatioTemporalIndices"))
    runs <- parLapply(cl,kFold, function(f)runwithout(data,pl,conf,map, fold=f,doDetailedRep = 0, doValidation = 1))
    stopCluster(cl)
  }
  flag = do.call(c, lapply(runs, function(f) f$opt$convergence))
  if(max(flag)>0){
    warning(paste0("Run ", which(flag>0), " did not converge"))
  }

  print("Done with fitting model with all training sets")
  indexVal = list(c(as.vector(runs[[1]]$obj$env$ADreportIndex()$validation)))

  if(nCores==1){
    repList <- lapply(runs,function(f){
      rep = sdreport(f$obj,bias.correct = TRUE,bias.correct.control = list(split =indexVal, sd = FALSE),
                     skip.delta.method = TRUE)
      rl = as.list(rep,"Est. (bias.correct)",report = TRUE)
      return(rl)
    })
  }else{
    cl <-makeCluster(ceiling(nCores/5),outfile = "")
    obj = lapply(runs,function(f)f$obj)
    clusterExport(cl,varlist=c("obj","indexVal"),envir=environment())
    clusterEvalQ(cl, library("spatioTemporalIndices"))
    repList <- parLapply(cl,obj,function(f){
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

      rep = sdreport(f,bias.correct = TRUE,bias.correct.control = list(split =indexVal, sd = FALSE),
                     skip.delta.method = TRUE)
      rl = as.list(rep,"Est. (bias.correct)",report = TRUE)
      return(rl)
    })
    stopCluster(cl)
  }


  predMatrix = matrix(0,length(kFold),run$data$numberOfLengthGroups )
  obsMatrix = predMatrix

  for(y in 1:length(run$conf$years)){
    for(l in 1:run$data$numberOfLengthGroups){
      obsMatrix[y,l] = mean(run$data$fishObsMatrix[kFold[[y]],l])
      predMatrix[y,l] = repList[[y]]$validation[l]/ length(kFold[[y]])
    }
  }

  return(list(obsMatrix = obsMatrix, predMatrix = predMatrix,runs = runs,repList = repList))
}




#' pointsInStratum
#' @param
#' @return
#' @export
#' @examples
#' @return A vector with index of points inside stratas providede in strataNumbers
pointsInStratum<-function(strataNumbers,stratas,locations){

  locationsSpObject=SpatialPoints(locations,proj4string=CRS("+proj=longlat"))
  toReturn = NULL;
  for(i in strataNumbers){
    g1=rgeos::readWKT(stratas$V2[i])
    proj4string(g1)=CRS("+proj=longlat")
    toReturn = c(toReturn,which(!is.na(over(locationsSpObject,g1))))
  }
  return(toReturn)
}


