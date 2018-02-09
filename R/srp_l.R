#' The (simpler) Smooth-Rough Partition linear regression model fitting
#'
#' This function performs same as \code{\link{srp.c}} except that constrained functional coefficient is estimated as a linear function.
#'
#' The estimation procedure of Smooth-Rough Partition model is described in "Regularised forecasting via smooth-rough partitioning of the regression coefficients", H. Maeng and P. Fryzlewicz (2018), preprint.
#'
#' @param x A matrix you wish to fit Smooth-Rough Partition model. The dimension of row is the number of variables.
#' @param y A vector you wish to use as a response variable in case of regressing \code{y} on \code{x}. If \code{y} is missing, the response variable is obtained in design matrix \code{x}.
#' @param maxq An integer specifying the maximum value of number of unconstrained parameters which the model can have. The default is 10.
#' @param plot If true, it gives the plot of estimated regression coefficients.
#' @return
#' \item{muhat}{The estimator of constant parameter.}
#' \item{bhat}{The vector of evaluated constrained (linear) functional regression coefficient.}
#' \item{ahat}{The vector of unconstrained regression coefficcient estimators.}
#' \item{etahat}{The vector containing both \code{bhat} and \code{ahat} with unevaluated form.}
#' \item{yhat}{The vector of estimated response variable.}
#' \item{SIC}{The vector of Schwarz criterion with length \code{maxq} which is computed for different number of unconstrained parameters.}
#' \item{qhat}{The optimal number of unconstrained parameters selected in the model.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{srp.c}}
#' @examples
#' x <- matrix(rnorm(10000), ncol=100)
#' srp.l(x)
#' @export

srp.l <- function(x, y, maxq=10, plot=T){

  ##### 1. Check the input
  if(!is.matrix(x))
    stop("The input 'x' should be a matrix.")
  if(missing(y)) {
    n <- dim(x)[2] # number of observation
    p <- dim(x)[1] # number of grid points each observation has
    y <- matrix(x[p, , drop=F], ncol=1) # the last column is picked for the response variable
    x <- x[c(1:p-1), , drop=F] # the predictors
  } else {
    if(!is.vector(y))
      stop("The input 'y' should be a vector.")
    y <- matrix(y, ncol=1) # the last column is picked for the response variable
    x <- x
    n <- dim(x)[2] # number of observation
    p <- dim(x)[1] + 1
  }

  ##### 2. Natural cubic spline interpolation of x
  cf0 <- ncs(x)$cf0

  ### define the basis for x and beta
  x.basis <- fda::as.fd(stats::splinefun((seq(0, 1, length.out=p-1)*(p-1)), x[, 1], method="natural"))$basis # equal basis for all observations
  M.basis <- fda::create.monomial.basis(rangeval=c(0, p-1), nbasis=2) # fit a polynomial(linear) basis

  #####################################################################################################
  ##### 3. Fitting the model SRP.L

  est_q <- sic.l(x.basis=x.basis, M.basis=M.basis, x=x, y=y, cf0=cf0, maxq=maxq, fixedq=F)
  qhat <- est_q$qhat
  sicq <- est_q$sicq

  fit <- sic.l(x.basis=x.basis, M.basis=M.basis, x=x, y=y, cf0=cf0, maxq=qhat, fixedq=T)

  if(plot==T){
    plot(c(1:(p-1)),c(fit$bhat, rep(NA, qhat)), type="l", lty=1, lwd=2
      ,ylim=range(c(fit$ahat, fit$bhat)), xlab=" ", ylab="",cex.axis=1.5)
    graphics::lines(c(1:(p-1)), c(rep(NA, p-1-qhat), fit$ahat),  type="b", lty=1, pch=19, lwd=2)
  }

  results <- list("muhat" = fit$muhat, "bhat" = fit$bhat, "ahat" = fit$ahat, "etahat" = fit$etahat, "yhat" = fit$yhat, "SIC" = sicq, "qhat" = qhat)
  return(results)


}
