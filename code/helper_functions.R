#' Collection of functions to assist in analyses
#' @author carlga ~ Carlos Gallardo
#'

#' Normalizes data column-wise using selected method
#' @param x A dataframe or matrix
#' @param len Vector with gene/transcript lengths (for RPKM/TPM normalization)
#' @param method Method used for normalization
#' @return Data frame with normalized data
#' 
#' @author carlga
#'
normalizeData <- function(x, len = NULL, method = NULL) {
  
  if(is.null(method)) stop('Please, indicate a valid method\n')
  
  x <- as.matrix(x)
  
  if(method == "CPM") {
    res <- t(t(x) * 1e6 / colSums(x))
  }
  else if(method == "RPKM") {
    if(length(len) != nrow(x)) stop('Length vector not matching gene number\n')
    res <- t(t(x) * 1e6 / colSums(x)) / (len / 1e3)
  }
  else if(method == "TPM") {
    if(length(len) != nrow(x)) stop('Length vector not matching gene number\n')
    x <- x / (len / 1e3)
    res <- t(t(x) * 1e6 / colSums(x))
  }
  else stop('This method is not available\n')
  
  return(as.data.frame(res))
}

#' Performs principal component analysis
#' @param x A dataframe or matrix
#' @param top_var Number of top variable features to use
#' @param center Center data
#' @param scale. Scale data
#' @return List with PCA results
#' 
#' @author carlga
#' 
doPCA <- function(x, top_var = 500, center = T, scale. = F) {
  require(matrixStats)
  
  rvars <- matrixStats::rowVars(as.matrix(x))
  select <- order(rvars, decreasing = TRUE)[seq_len(min(top_var, length(rvars)))]
  pca <- prcomp(t(x[select, ]), center = center, scale. = scale.)
  
  res <- list()
  res[['perc_var']] <- pca$sdev^2/sum(pca$sdev^2)*100
  res[['pcs']] <- as.data.frame(pca$x)
  res[['rotation']] <- as.data.frame(pca$rotation)
  
  return(res)
}
