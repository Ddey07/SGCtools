#' @title Get predictors of latent variables from observed data
#'
#' @description BLUPs of latent variables from observed data and estimated latent correlation matrix
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables. Binary variables must be coded as 0,1. Ordinal variables with l ordered categories must be coded as 0, 1, 2, ... , (l-1).
#' @param type A vector of length p containing type of variables in \code{X}, must be in the format "cont", "trunc", "ord", "bin". If not supplied, our function tries to guess the data-type from the observations.
#' @param lat.cov.est The estimate of the latent correlation matrix from \code{fromXtoRMixed}
#' @param impute.missing A logical variable indicating whether the missing data will be imputed with Missing at Random assumptions or not.
#' @return \code{getLatentPreds} returns a n by p matrix containing continuous predictions of observed mixed latent variables.
#' @references
#' Dey ., Zipunnikov V. (2022) "Semiparametric Gaussian Copula Regression modeling for Mixed Data Types (SGCRM)" <arXiv:2205.06868>
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/fromXtoR_ex.R
getLatentPreds = function(X, type = NULL, lat.cov.est, impute.missing = TRUE){
  PRx = preprocess_data(X)
  PR=PRx$features
  L0 = X
  for(i in 1:nrow(X)){
    L0[i,]=recover_row(i,X=X,lat.cov = lat.cov.est, feature=PR, type = type, impute = impute.missing)
  }
  return(L0)
}

