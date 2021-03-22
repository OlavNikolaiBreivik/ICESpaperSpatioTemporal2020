#' Initialie values for all model parameters and random effects.
#' @param data Data
#' @param conf Configurations
#' @details
#' @return a list containing initial values for all model parameters and random effects in the model.
#' @export
setPar <- function(data,conf){

  if(conf$lengthGroupsReduced[1] == conf$lengthGroupsReduced[2]){ #Latent effect when reducing length dimension
    xS = array(0.0, dim = c(dim(data$A_ListS[[1]])[2],NumberOfLengthGroups=max(conf$lengthGroupsReduced)+1))
    xST = array(0.0, dim = c(dim(data$A_ListST[[1]])[2],length(conf$years),NumberOfLengthGroups=max(conf$lengthGroupsReduced)+1))
  }else{
    xS = array(0.0, dim = c(dim(data$A_ListS[[1]])[2],NumberOfLengthGroups=max(conf$lengthGroupsReduced)))
    xST = array(0.0, dim = c(dim(data$A_ListST[[1]])[2],length(conf$years),NumberOfLengthGroups=max(conf$lengthGroupsReduced)))
  }

  nugget = array(0.0, dim = c(dim(data$fishObsMatrix)[2],dim(data$fishObsMatrix)[1]))


  nuggetIndex = array(0.0, dim = c(dim(data$fishObsMatrix)[2],length(data$xInt)))


  beta0 = array(0,dim=c(length(conf$years), length(conf$lengthGroups)))


  parameters <- list(beta0 = beta0,
                     betahaul =rep(1,length(conf$lengthGroups)),
                     betaSun =rep(0,max(1,conf$sunAlt[1])*4),
                     betaDepth = rep(0,data$Sdim*2),
                     log_lambda =c(5,5),
                     log_sigma =c(0,0,0),
                     log_kappa =c(-5,-5),
                     logSize =0.2,
                     tan_rho_t =1,
                     tan_rho_l =rep(1,3),
                     delta_z = rep(0,9),
                     xS = xS,
                     xST  = xST,
                     nugget = nugget,
                     nuggetIndex = nuggetIndex)
  return(parameters)
}


