#' The Smooth-Rough Partition model fitting
#'
#' This function performs the Smooth-Rough Partition linear regression with the input matrix.
#'
#' The estimation procedure of Smooth-Rough Partition model is described in "Regularised forecasting via smooth-rough partitioning of the regression coefficients", H. Maeng and P. Fryzlewicz (2018), preprint.
#'
#' @param x A matrix you wish to fit Smooth-Rough Partition model. The dimension of row is the number of variables which are pre-ordered in terms of their importance in prediction.
#' @param y A vector you wish to use as a response variable in case of regressing \code{y} on \code{x}. If \code{y} is missing, the response variable is obtained from the last row of \code{x}.
#' @param maxq An integer specifying the maximum number of unconstrained parameters which the model can have. The default is max(30, ceiling(0.1*dim(x)[1])).
#' @param L An integer specifying the dimension of b-spline expansion for the constrained (smoothed) parameters. The default is 35.
#' @param norder An integer specifying the order of b-splines. The default of 4 performs cubic splines.
#' @param inisp An initial value for optimising the tuning parameters and the default is 1.
#' @param plot If true, it gives the plot of estimated regression coefficients.
#' @return
#' \item{muhat}{The estimator of constant parameter.}
#' \item{bhat}{The vector of evaluated constrained functional regression coefficient.}
#' \item{ahat}{The vector of unconstrained regression coefficient estimators.}
#' \item{etahat}{The vector containing both \code{bhat} and \code{ahat} with unevaluated form.}
#' \item{yhat}{The vector of estimated response variable.}
#' \item{SIC}{The vector of Schwarz criterion with length \code{maxq} which is computed for the different number of unconstrained parameters.}
#' \item{qhat}{The optimal number of unconstrained parameters selected in the model.}
#' \item{sp}{The vector of two tuning parameters estimated by minimising generalised cross validation (GCV).}
#' \item{L}{The number of bases used for constrained regression parameters.}
#' \item{norder}{The order of b-splines specified.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{sic.c}}, \code{\link{predict.srp.c}}, \code{\link{srp.l}}
#' @examples
#' x <- matrix(rnorm(10000), ncol=100)
#' srp.c(x)
#' @export

srp.c <- function(x, y, maxq=max(30, ceiling(0.1*dim(x)[1])), L=35, norder=4, inisp=1, plot=T){

  ##### 1. Check the input
  if(!is.matrix(x))
    stop("The input 'x' should be a matrix.")
  if(maxq >= dim(x)[2]-L){
    L <- dim(x)[2] - maxq
  }
    #stop("'L' + 'maxq' should be less than the number of observation.")
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
  B.basis <- fda::create.bspline.basis(rangeval=c(0, p-1), norder=norder, nbasis=L)

  #####################################################################################################
  ##### 3. Fitting the model SRP.C

  est_q <- sic.c(x.basis=x.basis, B.basis=B.basis, x=x, y=y, cf0=cf0, maxq=maxq, fixedq=F, L=L, inisp=inisp)
  qhat <- est_q$qhat
  sicq <- est_q$sicq

  fit <- sic.c(x.basis=x.basis, B.basis=B.basis, x=x, y=y, cf0=cf0, maxq=qhat, fixedq=T, L=L, inisp=inisp)

  if(plot==T){
    plot(c(1:(p-1)),c(fit$bhat, rep(NA, qhat)), type="l", lty=1, lwd=2
      ,ylim=range(c(fit$ahat, fit$bhat)), xlab=" ", ylab="",cex.axis=1.5)
    graphics::lines(c(1:(p-1)), c(rep(NA, p-1-qhat), fit$ahat),  type="b", lty=1, pch=19, lwd=2)
  }

  results <- list("muhat" = fit$muhat, "bhat" = fit$bhat, "ahat" = fit$ahat, "etahat" = fit$etahat, "yhat" = fit$yhat, "SIC" = sicq, "qhat" = qhat, "sp" = fit$sp, "L" = fit$L, "norder" = norder)

  class(results) <- "srp.c"
  return(results)

}
