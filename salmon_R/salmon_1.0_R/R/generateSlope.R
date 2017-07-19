#' Slope Matrix
#' 
#' This function calculates a slope matrix from time-series log2FC expression data
#' 
#' @param lfc The numeric matrix or data.frame of log2FC data. Each row represents a gene and each column represents a sample.
#' @param tp A vector of time points of the samples in the matrix lfc. The length of the vector should be the same as the number of samples (i.e. the number of columns in the matrix \code{lfc}). 
#' @param group A vector of indices indicating the set of samples from a particular drug/compound treatment. The (time-series) samples from the same drug treatment experiment should have the same unique index. The length of the vector should be the same as the number of samples.  
#' 
#' @return \item{slopeMat}{the slope matrix}
#'
#' @details The time-series samples from the same experiment should have the same group index. A repitition experiment should have a different group index. Together with the information of time points \code{tp}, this information is needed to calculate time slope matrix.
#' 
#' For more than two time points in an experiment, the \code{2nd order accuracy finite difference} is used for calculating the slopes of the corresponding samples. If only two time points are available in an experiment, the \code{linear slope} from the two time points is used for that experiment.
#' 
#' 
#' @export
generateSlope <- function(lfc,tp,group){
  
  if(length(group)!=dim(lfc)[2]||length(tp)!=dim(lfc)[2]){
    stop("Dimension of group or timepoint does not match number of samples of the expression data.")
  }
  
  grp_ord <- order(group)
  group <- group[grp_ord]
  tp <- tp[grp_ord]
  lfc <- lfc[,grp_ord]
  
  grp_unique <- unique(group)
  if(length(grp_unique)==length(group)){
    stop("Slope matrix cannot be calculated because each group has only one sample.")
  }
  
  slopeMat <- vector(mode="numeric",length = 0L)
  #findiff <- dget("findiff.R")
  for (i in grp_unique) {
    grp_index <- which(group==i)
    if(length(unique(tp[grp_index]))!=1){
      if(length(unique(tp[grp_index]))!=length(tp[grp_index])){
        stop("Replicated samples of some time points in one group. Samples in each group shoud have distinct time points or the same time point.")
      }
      tp_order <- order(tp[grp_index])
      grp_index_order <- grp_index[tp_order]
      tp[grp_index] <- tp[grp_index_order]
      lfc[,grp_index] <- lfc[,grp_index_order]
        
      current_tp <- tp[grp_index]
      current_lfc <- lfc[,grp_index]
        
      l <- length(current_tp)
      if(l>2){
        for (j in 1:l) {
          if(j==1){
            slop_sp <- findiff(current_tp[1:3],current_lfc[,1:3],1)
          }else if(j==l){
            slop_sp <- findiff(current_tp[(l-2):l],current_lfc[,(l-2):l],3)
          }else{
            slop_sp <- findiff(current_tp[(j-1):(j+1)],current_lfc[,(j-1):(j+1)],2)
          }
          slopeMat <- cbind(slopeMat,slop_sp)
        }
      }else{
        slop_sp <- (current_lfc[,2]-current_lfc[,1])/(current_tp[2]-current_tp[1])
        slopeMat <- cbind(slopeMat,slop_sp)
      }
    }
  }
  return(slopeMat)
}