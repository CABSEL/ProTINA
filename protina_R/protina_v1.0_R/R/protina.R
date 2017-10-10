#' ProTINA
#' 
#' The main function for ProTINA for generating the protein scores for each drug treatment. 
#' @import glmnet BSDA foreach doParallel
#' 
#' @param lfc The numeric matrix or data.frame of log2FC data. Each row represents a gene and each column represents a sample.
#' @param slope The slope matrix from log2FC data. This matrix can be obtained using the function \code{generateSlope()}. If the data are not time-series, set slope to an empty matrix (i.e. \code{slope}=\code{NULL}). 
#' @param pgn The adjacency matrix of the protein-gene regulation network. This matrix can be created using the function \code{generatePGN()}. 
#' @param grplist The group index for protein scoring. This vector defines the samples for which the protein scores are computed. The length of this vector should be the same as the number of samples in the log2FC matrix. A single (aggregate) protein score is generated for samples with the same index. The group indices should be a consecutive integer starting from 1 to the number of groups.
#' @param kfold The number of folds used in the k-fold cross validation.
#' @param par A Boolean variable \code{TRUE} or \code{FALSE} indicating whether to use parallel computing. The default is \code{FALSE} (no parallel computation).
#' @param numCores The number of CPU cores to be used for parallel computing. This parameter is considered only if par is \code{TRUE}. The default is 4.
#' 
#' @return a list of two matrices
#' \item{Pscore}{The matrix of protein scores. Each row corresponds to a gene following the same order as the one in the log2FC data, while each column corresponds to a group of samples as defined in the grplist.}
#' \item{A}{The matrix of edge weights of the PGN. The (i,j)-th element of the matrix corresponds to the weight of the regulatory edge from protein j to gene i in the PGN, i.e.the regulation of the expression of gene i  by protein j. The rwos and columns of the matrix correspond to genes and proteins in the PGN.}
#' 
#' @export
protina <- function(lfc,slope=NULL,pgn,grplist,kfold=10,par=FALSE,numCores=4){
  
  n <- dim(lfc)[1]
  m <- dim(lfc)[2]
  
  norm_vec <- function(x) sqrt(sum(x^2))
  
  if(!is.null(slope)){
    ms <- dim(slope)[2]
    f_lfc <- apply(lfc, 1, norm_vec)*sqrt(m)
    f_slope <- apply(slope, 1, norm_vec)*sqrt(ms)
    rho <- f_lfc/f_slope
    rho[which(rho==Inf)] <- 0
    
    X0 <- cbind(lfc,slope)
    Im <- rbind(diag(m),matrix(0,nrow = ms,ncol = m))
  }else{
    X0 <- lfc
    Im <- diag(m)
  }
  
  dg <- unique(pgn$j)
  
  ###################################################################
  #DeltaNeTS ridge regression
  ###################################################################
  if(par){
    library(doParallel)
    library(foreach)
    cl <- makeCluster(numCores,outfile='')
    registerDoParallel(cl)
    
    A <- foreach(j=1:length(dg),.combine = rbind,.packages = "glmnet")%dopar%{
      cat(sprintf("Infering weights in PGN: (%6.3f)%%\n",(100*j/length(dg))))
      ppar <- pgn$i[pgn$j==dg[j]]
      
      if(!is.null(slope)){
        X <- X0
        X[,m+1:ms] <- rho[dg[j]]*X[,m+1:ms] #scaling slope variables
        y <- X[dg[j],]
      }else{
        X <- X0
        y <- X[dg[j],]
      }
      
      
      if(length(ppar)==1){
        Xin <- cbind(X[ppar,],Im)
      }else{
        Xin <- cbind(t(X[ppar,]),Im)
      }
    
      cvfit <- cv.glmnet(Xin, y, alpha = 0, intercept = FALSE, nfolds = kfold)
      coef.cvfit <- coef.cv.glmnet(cvfit,s='lambda.min')
      Aj <- vector(mode =  "integer",length = dim(Xin)[2])
      Aj[coef.cvfit@i] <- coef.cvfit@x

      r <- vector(length = n)
      r[ppar] <- Aj[1:(length(Aj)-m)]
      r
    }
  }else{
    A <- matrix(0,nrow = length(dg), ncol = n)
    for (j in 1:length(dg)) {
      cat(sprintf("Infering weights in PGN: (%6.3f)%%\n",(100*j/length(dg))))
      ppar <- pgn$i[pgn$j==dg[j]]
   
      
      if(!is.null(slope)){
        X <- X0
        X[,m+1:ms] <- rho[dg[j]]*X[,m+1:ms] #scaling slope variables
        y <- X[dg[j],]
      }else{
        X <- X0
        y <- X[dg[j],]
      }
      
      
      if(length(ppar)==1){
        Xin <- cbind(X[ppar,],Im)
      }else{
        Xin <- cbind(t(X[ppar,]),Im)
      }
      
      cvfit <- cv.glmnet(Xin, y, family='gaussian',alpha = 0, intercept = FALSE, nfolds = kfold)
      coef.cvfit <- coef.cv.glmnet(cvfit,s='lambda.min')
      Aj <- vector(mode =  "integer",length = dim(Xin)[2])
      Aj[coef.cvfit@i] <- coef.cvfit@x
      A[j,ppar] <- Aj[1:(length(Aj)-m)]
    }
  }
  
  ##########################################################
  #Z-test
  ##########################################################
  if(!is.null(slope)){
    res <- matrix(0,nrow=length(dg),ncol=m+ms)
    for (j in 1:length(dg)){
      X <- X0
      X[,m+1:ms] <- rho[dg[j]]*X[,m+1:ms] #scaling slope variables
      y <- X[dg[j],]
      res[j,] <- y-A[j,]%*%X
    }
  }else{
    Y <- X0[dg,]
    res <- Y-A%*%X0
  }

  
  ugrplist <- unique(grplist)
  muni <- length(ugrplist)
  cat(sprintf("Calculating Z-scores...\n"))
  
  pval <- matrix(1,length(dg),muni)
  z <- matrix(0,length(dg),muni)
  
  for (j in 1:length(dg)) {
    r <- res[j,]
    mr <- length(r)
    for (gri in 1:muni) {
      si <- c(1:length(grplist))[which(grplist==ugrplist[gri])]
      htest <- z.test(r[si],mu=0,sigma.x=sqrt(sum(r[setdiff(1:mr,si)]^2)/(mr-1-length(si))),conf.level=0.95,alternative='greater')
      pval[j,gri] <- htest$p.value
      z[j,gri] <- htest$statistic
    }
  }
  
  
  ####################################################
  #Combine z-scores
  ####################################################
  cat(sprintf('Combining z-scores for proteins...\n'))
  
  Aw <- A/apply(abs(A), 1, max)
  proti <- sort(unique(pgn$i))
  zp <- matrix(0,nrow = length(proti), ncol = muni)
  
  for (j in 1:length(proti)) {
    tgi <- match(unique(pgn$j[pgn$i==proti[j]]),dg)
    w <- Aw[tgi,proti[j]]
    if(length(tgi)>1) zp[j,] <- apply(z[tgi,]*w,2,sum)/norm_vec(w)
    else zp[j,] <- z[tgi,]*w/norm_vec(w)
  }
  
  Pscore <- matrix(0,nrow = n, ncol = muni)
  Pscore[proti,] <- zp 
  
  Af <- matrix(0,nrow = n, ncol = n)
  Af[dg,] = A
  
  cat(sprintf('The list for Pscore and matrix A is returned.\n'))
  return(list(Pscore=Pscore,A=Af))
}