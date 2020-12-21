

#' readRstoxFiles
#'
#' creates a list containing the data files from stox.
#'
#'@param conf Configurations, only years to include used.
#'
#'@export
readRstoxFiles<-function(conf){
  path = "inst/extdata/rstoxProjects/NEAcod/" # Will later be changed depending on rstox framework
  setJavaMemoryNumber=6e10 #Default value, maybe modified later

  files=list.files(path = path) #Later to be changed
  yearsAvailble = do.call(c,lapply(files, function(f) as.numeric(gsub("[^0-9]", "",  f))))
  yearIndex = which(yearsAvailble %in% conf$years)

  data= lapply(yearIndex,function(f) getBaseline(files[[f]],startProcess = 1, endProcess = 10,
                                                       parlist=list(RegroupLengthDist =list(LengthDist = "Process(StationLengthDist)", LengthInterval = 5))))
  names(data)=as.character(yearsAvailble[yearIndex])
  return(data)
}

#' setupData
#'
#' Set up data used by model
#'
#' @description
#' @param data A list of the data object returned by getBaseline, where each element of the list is a year
#' @param conf Configurations, see \code{\link{defConf}} for details.
#' @param confPred Prediction configurations
#' @return
#' @export
#' @examples
setupData = function(data = NULL,conf,confPred = NULL,useProj = FALSE,...){

  if(is.null(data) & !useProj){
    dataDir <- system.file("extdata", package = "spatioTemporalIndices")
    load(paste0(dataDir,"/cod.RData"))
    data = cod
  }else if(is.null(data)){
    data = readRstoxFiles(conf = conf)
  }

  #Set up structure used for the SPDE-procedure
  if(conf$cutoff[3]==0){
    meshS=createMesh(cutoff = conf$cutoff[1], cbound = conf$cbound[1])$mesh
    meshST=meshS
  }else{
    meshS=createMesh(cutoff = conf$cutoff[1], cbound = conf$cbound[1])$mesh
    meshST=createMesh(cutoff = conf$cutoff[2], cbound = conf$cbound[2])$mesh
  }

  spdeS = inla.spde2.matern(meshS, alpha=2)
  spdeST = inla.spde2.matern(meshST, alpha=2)
  spdeMatricesS = spdeS$param.inla[c("M0","M1","M2")]
  spdeMatricesST = spdeST$param.inla[c("M0","M1","M2")]


  fishObsList=vector("list", length(conf$lengthGroups))
  A_ListS=list(rep(1,length(conf$years)))
  A_ListST=list(rep(1,length(conf$years)))

  counter=1
  nStationsEachYear=rep(0,length(conf$years))
  distance=NULL;station=NULL;yearObs = NULL;locObs=NULL;locObsLatLon = NULL;timeInDay = NULL;depth=NULL;sunAlt = NULL;sunSetting = NULL;date = NULL;
  for(year in conf$years){

    d = dataInYear(year=year,conf=conf,data=data)
    dL=d[which(d$LengthGroup==conf$minLength),]
    loc=cbind(dL$UTMX,dL$UTMY)
    locLatLon=cbind(dL$x,dL$y)
    nStationsEachYear[counter]=dim(loc)[1]
    A_ListS[[counter]]=inla.spde.make.A(meshS,loc)
    A_ListST[[counter]]=inla.spde.make.A(meshST,loc)
    counter=counter+1

    #Store variables
    distance=c(distance,dL$dist)
    yearObs = c(yearObs,dL$year)
    station=c(station,dL$Station)
    locObs = rbind(locObs,loc)
    locObsLatLon = rbind(locObsLatLon,locLatLon)
    timeInDay = c(timeInDay,dL$timeInDay)
    sunAlt = c(sunAlt,dL$sunAlt)
    sunSetting = c(sunSetting,dL$sunSetting)
    date =  c(date,dL$date)
    depth=c(depth,dL$depth)
    LengthGroupCounter=1
    for(l in conf$lengthGroups){
      dL=d[which(d$LengthGroup==l),]
      fishObsList[[LengthGroupCounter]]=c(fishObsList[[LengthGroupCounter]],dL$fishObs)
      LengthGroupCounter=LengthGroupCounter+1
    }
  }
  #If depth is missing:
  missingDepth=which(is.na(depth))
  if(sum(missingDepth>0)>0){
    for(i in 1:length(missingDepth)){
      distMatrix=stats::dist(locObs)
      distances=spDistsN1(locObs,locObs[missingDepth[i],])
      minimumDistanceIndex=which(order(distances)==2)
      depth[missingDepth]=depth[minimumDistanceIndex]
    }
  }
  depth[depth<100]=100
  depth[depth>400]=400


  #Based on the list of fish-observations, we create a matrix
  fishObsMatrix=t(do.call(rbind, fishObsList))
  colnames(fishObsMatrix) = as.character(seq(min(conf$lengthGroups),max(conf$lengthGroups),by = conf$dLength))
  rownames(fishObsMatrix) = station

  predMatrix = fishObsMatrix*0 #If an element is 1, it is not included in the likelihood

  #Define the data and parameters given to TMB----------
  if(conf$sunAlt[1] <2){
    sunAltFormula = as.formula(paste( "fishObsMatrix ~ sin+cos"))
    sunAltFormulaIntRep = as.formula(paste( "sunAlt ~ sin+cos"))
  }else if(conf$sunAlt[1] ==2){
    sunAltFormula = as.formula(paste( "fishObsMatrix ~ sin+cos + sin2 + cos2"))
    sunAltFormulaIntRep = as.formula(paste( "sunAlt ~ sin+cos + sin2 + cos2"))
  }else if(conf$sunAlt[1] ==3){
    sunAltFormula = as.formula(paste( "fishObsMatrix ~ sin+cos + sin2 + cos2 + sin3 + cos3"))
    sunAltFormulaIntRep = as.formula(paste( "sunAlt ~ sin+cos + sin2 + cos2 + sin3 + cos3"))
  }else{
    stop("To many basis in Fourier approximation of sun heigth effect.")
  }

  #Constructing sun altitude as a variable from 0 to 1
  sunAltTrans = sunAlt*0
  for(i in 1:length(sunAlt)){
    dateThis = as.numeric(strsplit(date[i],";")[[1]])
    altMax = -999;altMin = 999
    for(t in seq(0,24,by = 0.05)){
      hour = floor(t)
      min = (t-hour)*60
      alt = altOfSun(min=min, hour=hour, day=dateThis[3], month=dateThis[2],lat=locObsLatLon[i,2], lon=locObsLatLon[i,1])$alt.of.sun
      if(altMax<alt){
        altMax = alt
      }
      if(altMin>alt){
        altMin = alt
      }
    }
    if(sunSetting[i]==0){
      sunAltTrans[i] = (sunAlt[i] - altMin)/(-altMin +altMax) /2
    }else{
      sunAltTrans[i] = 1-(sunAlt[i] - altMin)/(-altMin +altMax) /2
    }
  }
  sunAltTrans[which(sunAltTrans>1)]=1
  sunAltTrans[which(sunAltTrans<0)]=0

  gamSetup_sunAlt=mgcv::gam(sunAltFormula,
                               data=data.frame(haulID=rownames(fishObsMatrix),sin= sin(sunAltTrans *2*pi),
                                               cos=cos(sunAltTrans*2*pi),cos2=cos(2*sunAltTrans*2*pi),cos3=cos(3*sunAltTrans*2*pi),cos4=cos(4*sunAltTrans*2*pi),cos5=cos(5*sunAltTrans*2*pi),
                                               sin2=sin(2*sunAltTrans*2*pi),sin3=sin(3*sunAltTrans*2*pi),sin4=sin(4*sunAltTrans*2*pi),sin5=sin(5*sunAltTrans*2*pi),
                                               fishObsMatrix=fishObsMatrix[,1]),fit=FALSE)

  gamSetup_depth=mgcv::gam(fishObsMatrix~s(depth,bs="cs",k = conf$splineDepth[1]),
                              data=data.frame(haulID=rownames(fishObsMatrix),depth=depth,
                                              fishObsMatrix=fishObsMatrix[,1]),fit=FALSE)

  X_sunAlt = gamSetup_sunAlt$X[,-1]
  X_depth = gamSetup_depth$X[,-1]
  S_depth=as(gamSetup_depth$smooth[[1]]$S[[1]], "dgTMatrix")
  Sdim=nrow(S_depth)

  gamSetup_sunAltReport = gam(sunAltFormulaIntRep,
                              data=data.frame(sunAlt=seq(0,1,by = 0.01),sin=sin(seq(0,1,by = 0.01) *2*pi),sin2=sin(seq(0,1,by = 0.01)*2*2*pi),sin3=sin(seq(0,1,by = 0.01)*3*2*pi),
                                              cos = cos(seq(0,1,by = 0.01)*2*pi),cos2 = cos(seq(0,1,by = 0.01)*2*2*pi),cos3=cos(seq(0,1,by = 0.01)*3*2*pi)),fit=FALSE)
  X_sunAltReport = gamSetup_sunAltReport$X[,-1]

  #Maximum sun height is used in index calculation
  maxSunAlt = 0.5
  gamSetup_sunAltIntegrate = gam(sunAltFormulaIntRep,
                                    data=data.frame(sunAlt=rep(maxSunAlt,2),sin=sin(rep(maxSunAlt,2) *2*pi),
                                                    cos = cos(rep(maxSunAlt,2)*2*pi),cos2=cos(2*rep(maxSunAlt,2)*2*pi),cos3=cos(3*rep(maxSunAlt,2)*2*pi),
                                                    sin2=sin(2*rep(maxSunAlt,2)*2*pi),sin3=sin(3*rep(maxSunAlt,2)*2*pi)),fit=FALSE)
  X_sunAltIntegrate = gamSetup_sunAltIntegrate$X[,-1]

  X_depthReport = PredictMat(gamSetup_depth$smooth[[1]],data = data.frame(depth=as.numeric(seq(min(depth),max(depth),by = 1))))

  #Weighting used when smooting length dimension in latent effect
  weigthLength = conf$lengthGroupsReduced*0
  nCollaps = sum(conf$lengthGroupsReduced==1)
  counter = 1
  current = 1
  for(i in 1:length(weigthLength)){
    if(current==conf$lengthGroupsReduced[i]){
      current = current +1
      weigthLength[i] = 1
      counter = 1
    }else{
      weigthLength[i] = 1-counter/nCollaps
      counter = counter +1
    }
  }

  #Extract areas of each strata
  dataDir <- system.file("extdata", package = "spatioTemporalIndices")
  strata=read.delim(paste0(dataDir,"/vintertokt_barentshavny.txt"),head=F)
  stratas= lapply(strata$V2 , function(f){rgeos::readWKT(f)})
  stratasUTM= lapply(stratas, function(f){
    proj4string(f)=CRS("+proj=longlat +datum=WGS84");
    f<-spTransform(f,CRS("+proj=utm +zone=35 +datum=WGS84 +units=km +no_defs"));
    f})
  areas = rep(0,length(stratasUTM))
  for(s in 1:length(stratasUTM)){
    areas[s] = stratasUTM[[s]]@polygons[[1]]@area* (1/1.852)^2
  }

  #Include observations as a vector (needed for keep functionality)
  obsVector = apply(fishObsMatrix, 1, rbind)
  idxStart = seq(0,length(fishObsMatrix), by = dim(fishObsMatrix)[2])

  data <- list(A_ListS = A_ListS,
               A_ListST = A_ListST,
               dist = distance,
               fishObsMatrix = fishObsMatrix,
               obsVector = obsVector,
               idxStart = idxStart,
               numberOfLengthGroups=length(conf$lengthGroups),
               spdeMatricesS = spdeMatricesS,
               spdeMatricesST = spdeMatricesST,
               nStationsEachYear=nStationsEachYear,
               predMatrix = predMatrix,
               X_depth=X_depth,
               S_depth=S_depth,
               Sdim=Sdim,
               X_sunAlt = X_sunAlt,
               X_sunAltReport = X_sunAltReport,
               X_depthReport=X_depthReport,
               lengthGroupsReduced = conf$lengthGroupsReduced-1,
               weigthLength = weigthLength,
               areas = areas,
               nBasisSunAlt = max(1,conf$sunAlt[1]),
               obsModel = conf$obsModel,
               simulateProcedure = conf$simulateProcedure)

  attributes(data)$year = yearObs
  attributes(data)$meshS = meshS
  attributes(data)$meshST = meshST
  attributes(data)$locObs = locObs
  attributes(data)$locObsLatLon = locObsLatLon
  attributes(data)$depth = depth #Used when constructing integration points
  attributes(data)$X_sunAltIntegrate = X_sunAltIntegrate #Used when constructing integration points

  #Include integration points
  data = includeIntPoints(data,conf,confPred)

  data$selectedStratas = c(7,8,9,10,13,14,15,16,17,20) #Selected stratas in paper

  return(data)
}


#' dataInYear
#'
#' Prepares data from Rstox in a single year obtained by getBaseLine for the modeling routine.
#'
#' @description
#' @param year Year which is setup
#' @param conf Configurations
#' @param data Data from all years
#' @return data-object
#' @export
#' @examples
dataInYear = function(year,conf,data){

  g = data[[as.character(year)]]

  d = g$outputData$RegroupLengthDist #Use currently the grouped length dist. These are 5 cm length groups.
  d$WeightedCount[is.na(d$WeightedCount)] = 0

  for(id in unique(d$Station)){
    maxOnStation=max(d$LengthGroup[d$Station ==id])
    if(maxOnStation>conf$maxLength){
      #Create pluss group
      d$WeightedCount[d$Station ==id & d$LengthGroup == conf$maxLength] =
        sum(d$WeightedCount[!is.na(d$WeightedCount) & d$LengthGroup>=conf$maxLength &d$Station ==id ])

      #Remove observations above maxlength
      d = d[-which(d$Station ==id & d$LengthGroup > conf$maxLength),]
    }else if(maxOnStation==conf$maxLength){
      #Do nothing
    }else{
      #Add zeroes
      NumberOfNewLengthGroups=(conf$maxLength-maxOnStation) %/% 5
      newData=data.frame(SpecCat=rep(d$SpecCat[1],NumberOfNewLengthGroups),Station=rep(id,NumberOfNewLengthGroups),LengthGroup=seq(from=maxOnStation+5,to=conf$maxLength,by=5),LengthInterval=rep(5,NumberOfNewLengthGroups),WeightedCount=rep(0,NumberOfNewLengthGroups),LengthDistType=rep(d$LengthDistType[1],NumberOfNewLengthGroups))
      d=rbind(d,newData)
    }

    #Create minus group, probably never used
    minOnStation=min(d$LengthGroup[d$Station ==id])
    if(minOnStation<conf$minLength){
      d$WeightedCount[d$Station ==id & d$LengthGroup == conf$minLength] =
        sum(d$WeightedCount[!is.na(d$WeightedCount) & d$LengthGroup<= conf$minLength &d$Station ==id ])
      d = d[-which(d$Station ==id & d$LengthGroup < conf$minLength),]
    }else if(minOnStation==conf$minLength){
      #Do nothing
    }else{
      NumberOfNewLengthGroups=(minOnStation-conf$minLength) %/% 5
      newData=data.frame(SpecCat=rep(d$SpecCat[1],NumberOfNewLengthGroups),Station=rep(id,NumberOfNewLengthGroups),LengthGroup=seq(from=conf$minLength,to=minOnStation-5,by=5),LengthInterval=rep(5,NumberOfNewLengthGroups),WeightedCount=rep(0,NumberOfNewLengthGroups),LengthDistType=rep(d$LengthDistType[1],NumberOfNewLengthGroups))
      d=rbind(newData,d)
    }
  }

  #Merge data----------------

  #Define data colums
  nObs = length(d$SpecCat)
  d$x = rep(0,nObs)
  d$y = rep(0,nObs)
  d$cruise=rep(0,nObs)
  d$fishObs = rep(0,nObs)
  d$year=rep(year,nObs)
  d$dist = rep(0,nObs)
  d$starttime=rep(0,nObs)
  d$timeInDay=rep(0,nObs)
  d$sunAlt=rep(0,nObs)
  d$sunSetting=rep(0,nObs)
  d$date=rep(0,nObs)
  d$depth=rep(0,nObs)

  dS = g$outputData$FilterBiotic$FilterBiotic_BioticData_FishStation.txt
  d=d[order(d$Station,d$LengthGroup),]

  for(i in 1:nObs){
    #Extract the station number
    tmp = strsplit(d$Station[i],split = "/")[[1]][2]

    #Extract row of interest
    haul = match(tmp,dS$serialno)

    #Add data
    d$y[i] = dS$latitudestart[haul]
    d$x[i] = dS$longitudestart[haul]

    #Find sun height
    start = dS$starttime[haul]
    start = as.numeric(strsplit(start,":")[[1]])
    date = as.numeric(strsplit(dS$stopdate[haul],split = "/")[[1]])
    min = start[2];hour = start[1]; day = date[1]; month= date[2]; lat =  dS$latitudestart[haul]; lon =dS$longitudestart[haul]
    d$sunAlt[i] = altOfSun(min=min, hour=hour, day=day, month=month,lat=lat, lon=lon)$alt.of.sun
    d$date[i] = paste0(year,";",month,";",day)
    #Find if the sun is setting or not
    if(min>59.8){
      min = 59.8 #minor quick fix to not investigate next hour below.
    }
    if(altOfSun(min=min, hour=hour, day=day, month=month,lat=lat, lon=lon)$alt.of.sun >= altOfSun(min=min+0.1, hour=hour, day=day, month=month,lat=lat, lon=lon)$alt.of.sun){
      d$sunSetting[i] = 1
    }

    d$cruise[i]=dS$cruise[haul]
    d$fishObs[i] = round(d$WeightedCount[i]*dS$distance[haul])
    d$dist[i] = dS$distance[haul]

    #Include bottom depth: An average of maximum and minimum (if both exist)
    if(!is.na(dS$fishingdepthmax[haul]) &!is.na(dS$fishingdepthmin[haul])){
      d$depth[i] = (dS$fishingdepthmax[haul] +dS$fishingdepthmin[haul])/2
    }else if(!is.na(dS$fishingdepthmax[haul]) & is.na(dS$fishingdepthmin[haul])){
      d$depth[i]=dS$fishingdepthmax[haul]
    }else if(is.na(dS$fishingdepthmax[haul]) & !is.na(dS$fishingdepthmin[haul])){
      d$depth[i]=dS$fishingdepthmin[haul]
    }else{
      d$depth[i]=NA
    }
  }

  #NA observations correspond to 0 observations
  d$fishObs[is.na(d$fishObs)] = 0
  #---------------------------

  #Convert to UTM coordinates-------------------------------------------------------------------------
  loc = data.frame(d$x,d$y)
  names(loc) = c("X","Y")
  attr(loc, "projection") = "LL"
  attr(loc, "zone") = 35
  ddpcr::quiet(locUTM <- PBSmapping::convUL(loc))
  colnames(locUTM) = c("UTMX", "UTMY")

  d$UTMX = locUTM$UTMX
  d$UTMY = locUTM$UTMY
  #----------------------------------------------------------------------------------------------------
  return(d)
}





#' includeIntPoints
#' @return
#' @export
#' @examples
#' @return
includeIntPoints<-function(data,conf,confPred){
  points = constructIntPoints(confPred)
  ApredS = inla.spde.make.A(attributes(data)$meshS,loc = as.matrix(points$locUTM))
  ApredST = inla.spde.make.A(attributes(data)$meshST,loc = as.matrix(points$locUTM))

  data$ApredS = ApredS
  data$ApredST = ApredST
  data$predLoc = points$predLoc
  data$nStrata = length(unique(confPred$Strata))

  data$xInt = points$locUTM[,1]
  data$yInt = points$locUTM[,2]

  #Find depth covariate
  obs = SpatialPoints(attributes(data)$locObs,CRS("+proj=utm +zone=35 +datum=WGS84 +units=km +no_defs"))
  intPoints = SpatialPoints(points$locUTM,CRS("+proj=utm +zone=35 +datum=WGS84 +units=km +no_defs"))

  dist = gDistance(obs,intPoints, byid=T)
  minDist <- apply(dist, 1, function(x) order(x, decreasing=F)[2])

  xDepth = matrix(0,dim(points$locUTM)[1],dim(data$X_depth)[2])
  xDepth = data$X_depth[minDist,]
  data$X_depth_int = xDepth

  #Set sun altitude given at maximum
  X_sunAltIntegrateTmp = attributes(data)$X_sunAltIntegrate
  X_sunAltIntegrate = matrix(0,dim(points$locUTM)[1],dim(X_sunAltIntegrateTmp)[2])
  for(i in 1:dim(points$locUTM)[1]){
    X_sunAltIntegrate[i,1:dim(X_sunAltIntegrateTmp)[2]] = X_sunAltIntegrateTmp[1,1:dim(X_sunAltIntegrateTmp)[2]]
  }
  data$X_sunAltIntegrate = X_sunAltIntegrate
  return(data)
}
