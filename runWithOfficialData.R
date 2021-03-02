library(spatioTemporalIndices)
years = 1994:2019
lengthGroups = seq(20,100,by = 5)
conf = defConf(years = years, lengthGroups = lengthGroups)
load("dataNorwegian.RData")

par = setPar(dat,conf)
map = setMap(par,conf)


start_time <- Sys.time()
run = fitModel(data = dat,par = par,conf = conf,map = map)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


#Plot results
plotResults(what = "sunAlt",run = run)
plotResults(run,what = "depth")
partable(run)


#Apply bias-correct
library(spatioTemporalIndices)
setMKLthreads(1) #not profiting much by using more cores
start_time <- Sys.time()
indexIndices = c(as.vector(run$obj$env$ADreportIndex()$muReportFullExp))
x <- seq_along(indexIndices)
split = split(indexIndices, ceiling(x/10)) #Do bias correct in chunks
rep <- sdreport(run$obj,bias.correct = TRUE,bias.correct.control = list(split =split, sd = FALSE),
                skip.delta.method = TRUE, ignore.parm.uncertainty = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

rl = as.list(rep, "Est. (bias.correct)",report = TRUE)
rlOrig = as.list(run$rep, "Est",report = TRUE)
matplot(run$conf$years,rl$muReportFullExp/ rlOrig$muReportFullExp,
        main = "", xlab = "",ylab = "" )
mtext(text="Year",cex=1.7,side=1,line=2.4)
mtext(text="Multiplication factor",cex=1.7,side=2,line=2.4)
mtext("Bias correction factor",  cex=2.5, line=1)
abline(h = 0)



#Jit
set.seed(1234)
jj = jit(run,njit = 10,ncores = 10)

#OSA residuals
options(future.availableCores.methods = "mc.cores")
options(mc.cores = 20)
nHualsRemove = dim(run$data$fishObsMatrix)[1]-50
mm = max(run$data$fishObsMatrix[-(1:nHualsRemove),])
res <- oneStepPredict(run$obj, observation.name ="obsVector",
                      data.term.indicator="keep",discrete = TRUE, method = "oneStepGeneric",
                      conditional = 1:(nHualsRemove*length(run$data$lengthGroups)), parallel = TRUE,
                      range = c(0,round(10*mm)),
                      reverse = TRUE)


qqnorm(res$residual,main = "Quantile-quantile",
       cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
abline(0,1)

acf(res$residual,main = " ", ylab = " ",
    cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
title(main="Autocorrelation",cex.main = 2, ylab = "Autocorrelation", cex.lab = 1.5)



#Plot spatio-termporal illustration
mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
legend = FALSE
counter = 1
l = 60
for(y in seq(1994,2019,by = 1)){
  if(y ==2019)legend = TRUE
  addX = FALSE
  addY = FALSE

  if(counter%%mfrow[2]==1){
    addY = TRUE
  }
  counter = counter +1
  if(y >2015){
    addX = TRUE
  }
  plotResults(run,what = c("ST",y,l),range = c(-7,4),legend = legend)
}

mtext(text="Eastern direction (km)",cex=1.3,side=1,line=2,outer=TRUE)
mtext(text="Northern direction (km)",cex=1.3,side=2,line=2.4,outer=TRUE)
mtext("Spatio-temporal variation", outer=TRUE,  cex=2, line=0.5)






#Indices
index = mCPUESimultanious(run,rep,biasCorrect = TRUE)

cohortPlot = index$predIndexFull
div = max(cohortPlot)/11
cohortPlot = cohortPlot/div
plot(-1,1,xlab = "Year", ylab = "Length (cm)",
     xlim = c(1,length(run$conf$years)), ylim = c(1,length(run$conf$lengthGroups)),
     yaxt="n",
     xaxt = "n",
     main = "Indices at length",
     cex.main = 2,cex.lab = 1.7, cex = 1.6)
for(y in 1:dim(cohortPlot)[1]){
  points(rep(y,dim(cohortPlot)[2]),1:dim(cohortPlot)[2], cex = cohortPlot[y,])
}
axis(1,1:length(run$conf$years),run$conf$years)
axis(2,1:length(run$conf$lengthGroups),run$conf$lengthGroups)



#Indices in eatern part, Note that bias correction of "muReportSelectedExp" is needed here
index = mCPUESimultanious(run,repSelected,biasCorrect = TRUE)
ll = c(1,5,9,13)
for(l in ll){
  x11()
  cohort = index$predIndexSelected
  cohortL = index$predIndexSelectedL95
  cohortU = index$predIndexSelectedU95

  plot(run$conf$years,cohort[,l],type = 'l', ylim = c(0, max(cohortU[,l])),
       main = paste0("Index at length ", run$conf$lengthGroups[l], "cm in eastern zone"),
       xlab = "Prediction year",
       ylab = "Index at length",
       cex.main = 2,cex.lab = 1.7, cex = 1.6,cex.axis = 1.5)
  lines(run$conf$years,cohortL[,l],type = 'l',lty = 2)
  lines(run$conf$years,cohortU[,l],type = 'l',lty = 2)
}






#COG plot in space
CC = run$rep$cov
ii = which(names(run$rep$value)=="COGx" | names(run$rep$value)=="COGy")
S = CC[ii,ii]

rl = as.list(rep,"Est. (bias.correct)", report = TRUE)

indexHelper = matrix(1:length(rl$COGx),length(run$conf$years),length(run$conf$lengthGroups) )

for(length in c(20,40,60,80)){
  x11()
  par(oma=c(4,4,4,1.5),mar=c(0.5,0.5,0.0,0.0))
  plot(-10,-10,ylim = c(7700,8600),xlim = c(100,1300),
       cex.lab = 1.5,cex.axis = 1.8, cex.main=2, cex.sub= 1.8,
       main = "",
       xlab = "",
       ylab = "")

  mtext(text="East (km)",cex=1.3,side=1,line=2,outer=TRUE)
  mtext(text="North (km)",cex=1.3,side=2,line=2.4,outer=TRUE)
  mtext(paste0("COG, length ", length, "cm") , outer=TRUE,  cex=2, line=0.5)

  newmap <- map("world", c("Norway","Sweden","Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
  mapTmp = data.frame(newmap$x,newmap$y)
  mapTmp[which(is.na(mapTmp[,1])),] = 3.141592 #Need something different from NA
  names(mapTmp) = c("X","Y")
  attr(mapTmp, "projection") = "LL"
  attr(mapTmp, "zone") = 35
  ddpcr::quiet(mapTmp <- PBSmapping::convUL(mapTmp))
  colnames(mapTmp) = c("UTMX", "UTMY")
  mapTmp[which(is.na(newmap$x)),] = NA
  polygon(mapTmp,col = 'lightgrey')

  col = rev(heat.colors(length(run$conf$years) + 10))

  l = which(run$conf$lengthGroups==length)

  for(y in seq(1,length(run$conf$years),by = 2)){
    tmp =indexHelper[y,l]
    mu = c(rl$COGx[tmp], rl$COGy[tmp])
    testS = matrix(NA,ncol = 2, nrow= 2)
    testS[1,1] = S[tmp,tmp]
    testS[2,2] = S[length(rl$COGx) +tmp+1,length(rl$COGx) +tmp+1]

    testS[1,2] = S[tmp,length(rl$COGx) +tmp+1]
    testS[2,1] = S[length(rl$COGx) +tmp+1,tmp]
    lines( ellipse::ellipse( testS, centre = mu,level = 0.95) , col=col[y+10],type = 'l',lw = 2)
    Sys.sleep(0.01)
  }

  if(length==80){
    x<- rep(1,10)#Dummy variable to get aruments on legend
    y<- 1994:2016
    z<- outer( x,y,"+")

    image.plot( x,y,z,legend.only=TRUE, col = col[10:length(col)],
                legend.args=list( text="Year",cex=2),
                smallplot = c(0.8, 0.85, .1, .9),
                axis.args=list(cex.axis=1.3))
  }
}

#Run with different latent effects---------------------
conf = run$conf

conf$spatial = 1
conf$spatioTemporal = 0
conf$nugget = 1
runS = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 1)

conf = run$conf
conf$spatial = 0
conf$spatioTemporal = 1
conf$nugget = 1
runST = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 1)

conf = run$conf
conf$spatial = 0
conf$spatioTemporal = 0
conf$nugget = 1
runNoSST = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 1)

conf = run$conf
conf$spatial = 1
conf$spatioTemporal = 0
conf$nugget = 0
runOnlyS = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 1)

conf = run$conf
conf$spatial = 0
conf$spatioTemporal = 0
conf$nugget = 0
runNo = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 1)

AIC(runNo,runOnlyS,runNoSST,runS,runST,run)
#------------------------------------------------------------------------



#Different combinations of parameters
conf = run$conf
conf$splineDepth = c(6,1)
conf$sunAlt = c(0,0)
runD = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 2)

conf$splineDepth = c(6,0)
conf$sunAlt = c(1,0)
runS = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 2)

conf$splineDepth = c(6,0)
conf$sunAlt = c(0,0)
runN = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 2)

AIC(runD,runS,runN,run)




#Length dependent covariates
conf = run$conf
conf$splineDepth = c(6,2)
conf$sunAlt = c(1,1)
runLD = fitModel(data = run$data,par = run$par,conf = conf,doDetailedRep = 2)

plotResults(what = "sunAlt",run = runLD)
plotResults(runLD,what = "depth")




#Run without storing sdrep, memory issues in validation when applying several cores.
runVal = fitModel(data = dat,par = par,conf = conf,map = map, doDetailedRep = 0)

start_time <- Sys.time()
val = validateModel(runVal,nCores=8)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


plot(log(val$obsMatrix),log(val$predmatrix),
     main = "Predictions vs observations",xlab = "Observation",ylab = "Prediction",
     cex.lab = 1.5, cex.main = 1.7,cex.axis = 1.4,
     ylim = c(-4.5,4.5))
abline(0,1,lw = 3)
mean(abs(log(valSmall$predmatrix)-log(valSmall$obsMatrix)))
mean((log(valSmall$predmatrix)-log(valSmall$obsMatrix))^2)
