

#' effetctiveFishingWidth (not used in paper)
#'
#' Returns the effective fishing width for a specific 5 cm length group.
#' @param length in cm
#' @return
#' @export
#' @examples
#' @return Returns the length dependent catchability constant
effetctiveFishingWidth <- function(length){
  if(length<=15){
    a=5.91*15^0.43
  }else if(length<=60){
    a=5.91*length^0.43
  }else{
    a=5.91*62^0.43
  }

  scale = 1852/a
  return(scale)
}


#' adjustIndexLength (not used in paper)
#'
#' @param
#' @return
#' @export
#' @examples
#' @return
adjustIndexLength <- function(index,length = NULL,stratas = stratas){
  if(is.matrix(index)){
    for(l in 1:length(length)){
      index[,l] = index[,l]*effetctiveFishingWidth(length[l])
    }
    rownames(index) = stratas; colnames(index) = length
    return(index)
  }else if(is.vector(index)){
    for(l in 1:length(length)){
      index[l] = index[l]*effetctiveFishingWidth(length[l])
    }
    return(index)
  }else if(is.array(index)){
    for(l in 1:length(length)){
      index[,,l] = index[,,l]*effetctiveFishingWidth(length[l])
    }
    return(index)
  }
}




#' constructIntPoints
#' @return
#' @export
#' @examples
#' @return
constructIntPoints<-function(confPred){
  dataDir <- system.file("extdata", package = "spatioTemporalIndices")

  #Read stratas and convert to UTM
  preStratas=read.delim(paste0(dataDir,"/vintertokt_barentshavny.txt"),head=F)
  stratasPolygons = lapply(confPred$Strata,function(f) rgeos::readWKT(preStratas$V2[f]))
  for(i in 1:length(stratasPolygons)){
    stratasPolygons[[i]]@proj4string = CRS("+proj=longlat")
  }
  utmCRS = CRS("+proj=utm +zone=35 +datum=WGS84 +units=km +no_defs")
  stratasPolygonsUTM = lapply(stratasPolygons,function(f) spTransform(f,utmCRS))
  stratasPolygonsUTMtmp = lapply(1:length(stratasPolygonsUTM), function(x) {
    stratasPolygonsUTM[[x]]@polygons[[1]]@ID = as.character(x)
    stratasPolygonsUTM[[x]]
  })
  unionStratas = SpatialPolygons(lapply(stratasPolygonsUTMtmp, function(x){x@polygons[[1]]}))
  unionStratas@proj4string = utmCRS

  #Construct a grid of integration points
  points = sp::makegrid(unionStratas,confPred$nIntPoints)
  pointsSP = sp::SpatialPoints(points,utmCRS)
  inside = rep(1,dim(pointsSP@coords)[1])
  inside[which(is.na(over(pointsSP,unionStratas)))] = 0
  points = points[which(inside==1),]

  #Remove a few integration points outside of stratas when converting back to lat lon
  pointsSP = sp::SpatialPoints(points,utmCRS)
  pointsXY =  spTransform(pointsSP,CRS("+proj=longlat"))
  inside = rep(0,dim(pointsXY@coords)[1])
  joinedXY = spTransform(unionStratas,CRS("+proj=longlat"))
  inside[which(!is.na(over(pointsXY,joinedXY)))] = 1
  remove = which(inside==0)
  if(length(remove)>0){
      points = points[-remove,]
  }
  pointsSP = sp::SpatialPoints(points,utmCRS)

  #Define data frame with integration points to be returned
  locUTM = data.frame(points[,1],points[,2]) #To be returned
  colnames(locUTM) = c("UTMX", "UTMY")

  #Define matrix to link stratas with integration points
  predLoc = matrix(-1,length(confPred$Strata),2000)
  for(i in 1:length(confPred$Strata)){
    insideThis =  which(!is.na(over(pointsSP,stratasPolygonsUTM[[i]])))
    predLoc[i,1:length(insideThis)] = insideThis-1
  }


  return(list(locUTM = locUTM, predLoc = predLoc))
}


