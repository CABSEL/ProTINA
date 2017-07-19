#' Protein-Gene Network Construction
#' 
#'  Generate protein-gene regulatory network from gene regulatory network and protein protein interaction network
#' @import slam
#'  
#' @param glist The vector of genes in the same order of the genes in log2FC data. The length of GList should be the same as the number of rows in log2FC data.
#' 
#' @param tftg The matrix of TF-gene interactions. The first column is the list of TFs, and the second column is the list of genes regulated by the corresponding TFs. The third column is optional, and if present, the column should contain the (confidence) score for each interaction.
#' 
#' @param ppi The matrix of protein-protein interactions. Each row of the first two columns give the protein pairs with interactions. The third column is optional, and if present, the column should contain the (confidence) score for each interaction.
#' 
#' @param tftg_thre A threshold for TF-gene interactions. This variable is used only when the confidence score of TF-gene interactions are given in the matrix tftg. Any TF-gene interactions with confidence scores lower than the threshold will be excluded.
#' 
#' @param ptf_thre A threshold for protein-TF interaction. This variable is used only when the confidence score of protein-TF interactions are given in the matrix ppi.  Any protein-TF interactions with the scores lower than the threshold will be excluded.
#' 
#' @param ppi_thre A threshold for protein-protein interaction. This variable is used only when the confidence score of protein-protein interactions are given in the matrix ppi.  Any protein-protein interactions with the scores lower than the threshold will be excluded.
#' 
#' @return \item{pgn}{The adjacency matrix of PGN}
#' 
#' @export
generatePGN <- function(glist,tftg,ppi,tftg_thre=NULL,ptf_thre=NULL,ppi_thre=NULL){
  glist <- sapply(glist, toupper)
  
  if (ncol(tftg)==3){
  names(tftg) <- c("tf","tg","score")
  }else names(tftg) <- c("tf","tg")
  if (ncol(ppi)==3){
    names(ppi) <- c("protein1","protein2","score")
  }else names(ppi) <- c("protein1","protein2")
  
  
  #GRN construction
  if(ncol(tftg)==3&&!is.null(tftg_thre)){
    tftg <- tftg[which(tftg$score>=tftg_thre),]
  }
  tftg <- as.data.frame(sapply(tftg, toupper))
  
  edge <- matrix(nrow = nrow(tftg),ncol = 3)
  edge[,1] <- match(tftg$tf,glist)
  edge[,2] <- match(tftg$tg,glist)
  if(ncol(tftg)==3){
    edge[,3] <- as.numeric(as.character(tftg$score))
  }else edge[,3] <- 1
  edge <- edge[!is.na(edge[,1])&!is.na(edge[,2]),]
  
  grn <- matrix(0,nrow = length(glist),ncol = length(glist))
  grn[edge[,1:2]] <- edge[,3]
  grn <- as.simple_triplet_matrix(grn)
  
  tfn <- length(unique(edge[,1]))
  dgn <- length(unique(edge[,2]))
  edgen <- length(grn$i)
  cat(sprintf("TF-Gene network with %d TFs, %d genes and %d interactions has been generated.\n",tfn,dgn,edgen))
  
  #PPI construction
  ppi <- as.data.frame(sapply(ppi, toupper))
  
  edge <- matrix(nrow = nrow(ppi),ncol = 3)
  edge[,1] <- match(ppi$protein1,glist)
  edge[,2] <- match(ppi$protein2,glist)
  if(ncol(ppi)==3){
    edge[,3] <- as.numeric(as.character(ppi$score))
  }else edge[,3] <- 1
  edge <- edge[which(edge[,1]!=edge[,2]&!is.na(edge[,1])&!is.na(edge[,2])),]
  
  first2columns <- as.data.frame(t(apply(edge[,1:2], 1, sort)))
  first2columns$V3 <- edge[,3]
  first2columns <- first2columns[order(-first2columns$V1,-first2columns$V2,-first2columns$V3),]
  dup <- duplicated(first2columns[,1:2])
  ppi_idx <- first2columns[!dup,]
  colnames(ppi_idx) <- c("V1","V2","V3")
  
  #Protein1-TF network
  tf_idx <- unique(grn$i)
  e1 <- ppi_idx[,1]%in%tf_idx|ppi_idx[,2]%in%tf_idx
  ppi_ly1 <- ppi_idx[e1,]
  proti1 <- setdiff(unique(c(ppi_ly1[,1],ppi_ly1[,2])),tf_idx)
  
  #set the order (prot1-TF)
  x <- ppi_ly1[,1]%in%tf_idx
  ppi_ly1[x,1:2] <- ppi_ly1[x,2:1]
  
  #copy TF2-TF1 interaction
  copyi <- ppi_ly1[,1]%in%tf_idx&ppi_ly1[,2]%in%tf_idx
  copy <- ppi_ly1[copyi,c(2,1,3)]
  colnames(copy) <- c("V1","V2","V3")
  ppi_ly1 <- rbind(ppi_ly1,copy)
  if(!is.null(ptf_thre)){
    ppi_ly1 <- ppi_ly1[which(ppi_ly1[,3]>=ptf_thre),]
  }
  
  ppn1 <- matrix(0,nrow = length(glist),ncol = length(glist))
  ppn1[as.matrix(ppi_ly1[,1:2])] <- ppi_ly1[,3]
  ppn1 <- as.simple_triplet_matrix(ppn1)
  
  #summary
  prot1n <- length(unique(ppi_ly1[,1]))
  tfn <- length(unique(ppi_ly1[,2]))
  edgen <- length(ppn1$i)
  cat(sprintf("Protein-TF network with %d proteins, %d TFs and %d interactions has been generated.\n",prot1n,tfn,edgen))
  
  #Prot2-Prot1 network
  e2 <- ppi_idx[,1]%in%proti1|ppi_idx[,2]%in%proti1
  e2 <- setdiff(which(e2),which(e1))
  ppi_ly2 <- ppi_idx[e2,]
  proti2 <- setdiff(unique(c(ppi_ly2[,1],ppi_ly2[,2])),proti1)
  
  #set the order(Prot2-Prot1)
  y <- ppi_ly2[,1]%in%proti1
  ppi_ly2[y,1:2] <- ppi_ly2[y,2:1]
  
  #copy Prot1b-Prot1a interaction
  copyi <- ppi_ly2[,1]%in%proti1&ppi_ly2[,2]%in%proti1
  copy <- ppi_ly2[copyi,c(2,1,3)]
  colnames(copy) <- c("V1","V2","V3")
  ppi_ly2 <- rbind(ppi_ly2,copy)
  if(!is.null(ppi_thre)){
    ppi_ly2 <- ppi_ly2[which(ppi_ly2$V3>=ppi_thre),]
  }
  
  ppn2 <- matrix(0,nrow = length(glist),ncol = length(glist))
  ppn2[as.matrix(ppi_ly2[,1:2])] <- ppi_ly2[,3]
  ppn2 <- as.simple_triplet_matrix(ppn2)
  
  #summary
  prot2n <- length(unique(ppi_ly2[,1]))
  prot1n <- length(unique(ppi_ly2[,2]))
  edgen <- length(ppn2$i)
  cat(sprintf("Protein-protein network with %d upstream proteins, %d proteins and %d interactions has been generated.\n",prot2n,prot1n,edgen))
  
  #Protein-gene network construction
  pgn <- crossprod_simple_triplet_matrix(as.simple_triplet_matrix(t(ppn2+diag(length(glist)))),ppn1)
  pgn <- crossprod_simple_triplet_matrix(as.simple_triplet_matrix(t(pgn+diag(length(glist)))),grn)
  pgn <- as.simple_triplet_matrix(pgn)
  
  #summary
  protn <- length(unique(pgn$i))
  dgn <- length(unique(pgn$j))
  edgen <- length(pgn$v)
  cat(sprintf("Protein-gene network with %d proteins and TFs, %d downstream genes and %d interactions has been generated.\n",protn,dgn,edgen))
  
  return(pgn)
}