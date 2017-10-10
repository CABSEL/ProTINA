#' 2nd order accuracy finite difference from three points
#' 
#'  This function implements a 2nd order accurate finite difference for calculating slopes (time-derivatives) of the log2FC data using three time points. The function is used in \code{generateSlope} function.
#'  
#' @param tpoints Three consecutive time points in increasing order in an experiment
#' @param f Three columns of the log2FC matrix. Columns correspond to the samples given for \code{tpoints}.
#' @param tsi The indication of the current time point (1 for the initial time point, 2 for middle time points, and 3 for the last time point) 
#' @return \item{dfdt}{A vector of slopes at the current time point}
#' 
#' @details For the initial time point (\code{tsi}=1) in an experiment, \code{tpoints} are the set of the initial and next two time ptoins. For the last time point in the experiment (\code{tsi}=3),\code{tpoints} are the set of the last and previous two time ptoins in increasing order. For a middle time point (\code{tsi}=2), \code{tpoints} are a current time point and two neighbor time points in increasing order. 
#' Depending on \code{tsi}, three different strategies (forward, centeral, and backward) are used in approximation: forward for \code{tsi}=1, central for \code{tsi}=2, and backward for \code{tsi}=3.
#'  


findiff <- function(tpoints,f,tsi){
  tdiff <- tpoints - tpoints[tsi]
  if(length(which(tdiff>0))==0){
    deltsind <- max(which(tdiff<0))
    signdelts <- -1
  }else{
    deltsind <- min(which(tdiff>0))
    signdelts <- 1
  }
  delt <- tdiff/abs(tdiff[deltsind])

  Parind <- 1:length(tpoints)
  Parind <- Parind[-c(tsi,deltsind)]

  xmat <- delt[Parind]^2
  sol <- -1/xmat

  dfdt <- (-(1+sol)*f[,tsi]+f[,deltsind]+f[,Parind]*sol)/((signdelts+delt[Parind]*sol)*abs(tdiff[deltsind]))
  return(dfdt)
}
