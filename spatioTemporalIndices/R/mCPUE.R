

#' mCPUE
#' @param run Estimated model
#' @return
#' @export
#' @examples
#' @return Predicted index in all stratas
mCPUEstrata<-function(run){

  rl = as.list(run$rep,"Est",report = TRUE)
  rlSd = as.list(run$rep,"Std",report = TRUE)

  index = rl$muReport
  indexL95= rl$muReport - 1.96*rlSd$muReport
  indexU95= rl$muReport + 1.96*rlSd$muReport

  predIndex = exp(index)
  predIndexL95 = exp(indexL95)
  predIndexU95 = exp(indexU95)

  predIndexLAdj = adjustIndexLength(predIndex,length = run$conf$lengthGroups,stratas = run$confPred$Strata)
  predIndexLAdjL95 = adjustIndexLength(predIndexL95,length = run$conf$lengthGroups,stratas = run$confPred$Strata)
  predIndexLAdjU95 = adjustIndexLength(predIndexU95,length = run$conf$lengthGroups,stratas = run$confPred$Strata)

  return(list(predIndexLAdj = predIndexLAdj,predIndexLAdjL95=predIndexLAdjL95,predIndexLAdjU95 = predIndexLAdjU95  ))

}




#' mCPUESimultanious
#' @param run Estimated model
#' @return
#' @export
#' @examples
#' @return Estimated index in all stratas simulaniously, and in selected areas provided to the C-side
mCPUESimultanious<-function(run,rep = NULL, biasCorrect = FALSE){

  nYears = length(run$conf$years)
  nStrata = 26

  rl = as.list(run$rep,"Est",report = TRUE)
  rlSd = as.list(run$rep,"Std",report = TRUE)

  if(biasCorrect){
    if(is.null(rep)) stop("Need to provide sdrep with bias.correct = TRUE")
    rl = as.list(rep,"Est. (bias.correct)", report = TRUE)
  }

  indexFull = log(rl$muReportFullExp)
  sdFull = rlSd$muReportFull
  indexFullL95 = indexFull - 1.96*sdFull
  indexFullU95 = indexFull + 1.96*sdFull

  indexSelected = log(rl$muReportSelectedExp)
  sdSelected = rlSd$muReportSelected
  indexSelectedL95 = indexSelected - 1.96*sdSelected
  indexSelectedU95 = indexSelected + 1.96*sdSelected

  predIndexFull = exp(indexFull)
  predIndexFullL95 =exp(indexFullL95)
  predIndexFullU95 = exp(indexFullU95)

  predIndexSelected = exp(indexSelected)
  predIndexSelectedL95 = exp(indexSelectedL95)
  predIndexSelectedU95 = exp(indexSelectedU95)


  rownames(predIndexFull) = run$conf$years; colnames(predIndexFull) = run$conf$lengthGroups
  rownames(predIndexFullL95) = run$conf$years; colnames(predIndexFullL95) = run$conf$lengthGroups
  rownames(predIndexFullU95) = run$conf$years; colnames(predIndexFullU95) = run$conf$lengthGroups
  rownames(predIndexSelected) = run$conf$years; colnames(predIndexSelected) = run$conf$lengthGroups
  rownames(predIndexSelectedL95) = run$conf$years; colnames(predIndexSelectedL95) = run$conf$lengthGroups
  rownames(predIndexSelectedU95) = run$conf$years; colnames(predIndexSelectedU95) = run$conf$lengthGroups

  #Apply length dependent catchability (not used in paper)
  predIndexFullLAdj = adjustIndexLength(predIndexFull,length = run$conf$lengthGroups,stratas = run$confPred$Strata)
  predIndexFullLAdjL95 = adjustIndexLength(predIndexFullL95,length = run$conf$lengthGroups,stratas = run$confPred$Strata)
  predIndexFullLAdjU95 = adjustIndexLength(predIndexFullU95,length = run$conf$lengthGroups,stratas = run$confPred$Strata)

  predIndexSelectedLAdj = adjustIndexLength(predIndexSelected,length = run$conf$lengthGroups,stratas = run$confPred$Strata)
  predIndexSelectedLAdjL95 = adjustIndexLength(predIndexSelectedL95,length = run$conf$lengthGroups,stratas = run$confPred$Strata)
  predIndexSelectedLAdjU95 = adjustIndexLength(predIndexSelectedU95,length = run$conf$lengthGroups,stratas = run$confPred$Strata)

  return(list(indexFull = indexFull,indexSelected= indexSelected,
              sdFull = sdFull,sdSelected = sdSelected,
              predIndexFull = predIndexFull,predIndexFullL95 = predIndexFullL95,predIndexFullU95 = predIndexFullU95,
              predIndexSelected = predIndexSelected,predIndexSelectedL95 = predIndexSelectedL95,predIndexSelectedU95 = predIndexSelectedU95,
              predIndexFullLAdj = predIndexFullLAdj,predIndexFullLAdjL95 = predIndexFullLAdjL95,predIndexFullLAdjU95 = predIndexFullLAdjU95,
              predIndexSelectedLAdj = predIndexSelectedLAdj,predIndexSelectedLAdjL95 = predIndexSelectedLAdjL95,predIndexSelectedLAdjU95 = predIndexSelectedLAdjU95))

}



