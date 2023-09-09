#' @title Fast calculation of Kendall's Tau from observed data
#'
#' @description Calculate Kendall's Tau using fast approach of counting swaps by Samuel Perrault. Updated to handle missing data. If missing, Kendall's Tau is calculated based on pairwise complete observations.
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables.
#' @return \code{Kendall_mixed} returns a p by p matrix containing the Kendall's Tau of pairwise variables.
#' @references
#' Perreault S., "Efficient inference for Kendall's tau", <arXiv:2206.04019>.
#' @export
#' @example man/examples/fromXtoR_ex.R

Kendall_mixed = function(X){
  CDN <- conc_ties(X)
  C <- CDN$C
  D <-  CDN$D
  N <- CDN$N
  C_sum <- apply(C, c(2,3), function(x){sum(x, na.rm = TRUE)})
  D_sum <- apply(D, c(2,3), function(x){sum(x, na.rm = TRUE)})
  T2 <- (C_sum - D_sum)/(2*choose(N,2))
  return(T2)
}
