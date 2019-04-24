#' The natural cubic spline interpolation
#'
#' This function performs the natural cubic spline interpolation of a design matrix.
#'
#' Usually only called by \code{\link{srp.c}} and \code{\link{srp.l}}.
#'
#' @param x The design matrix.
#' @return
#' \item{cf0}{The coefficient matrix for B-splines obtained by natural cubic spline interpolation of \code{x}.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{srp.c}}, \code{\link{srp.l}}
#' @examples
#' x <- matrix(rnorm(100), ncol=10)
#' ncs(x)$cf0
#' @export

ncs <- function(x){
  ##### Natural cubic spline interpolation of x
  t <- seq(0, 1, length.out=dim(x)[1])*(dim(x)[1]) # the time interval in which the predictors live
  cf0 <- matrix(NA, nrow=dim(x)[1]+2, ncol=dim(x)[2]) # the coefficient matrix for the interpolated x
  ### interpolation of each observation x_i
  for(i in 1:(dim(x)[2])){
    cf0[, i] <- fda::as.fd(stats::splinefun(t, x[, i], method="natural"))$coefs
  }
  return(list("cf0"=cf0, "t"=t))
}
