#' The smoothness change-point detection of regression coefficients in linear model
#'
#' This function performs the optimisation of the number of unconstrained regression parameters in Smooth-Rough Partition model and gives the change-point of smoothness in regression parameters.
#'
#' The Smooth-Rough Partition model is described in "Regularised forecasting via smooth-rough partitioning of the regression coefficients", H. Maeng and P. Fryzlewicz (2018), preprint.
#'
#' @param x.basis The b-spline basis defined for interpolated x in \code{\link{srp.c}}
#' @param B.basis The b-spline basis defined for constrained regression coefficient.
#' @param x The design matrix used in \code{\link{srp.c}}.
#' @param y The response variable used in \code{\link{srp.c}}.
#' @param cf0 The coefficient matrix obtained by natural cubic spline interpolation of x in \code{\link{ncs}}.
#' @param maxq The maximum value of possible number of unconstrained parameters if \code{fixedq} is FALSE. Otherwise, it is considered as an unique number of unconstrained parameters.
#' @param fixedq If TRUE, \code{maxq} is considered as a fixed number of unconstrained parameters and if FALSE, \code{maxq} is a maximum and a sequence of possible values are investigated to select the optimal.
#' @param L The dimension of b-spline expansion for constrained parameters used in \code{\link{srp.c}}.
#' @return The following components are obtained only when \code{fixedq} is FALSE:
#' \item{qhat}{The optimal number of unconstrained parameters.}
#' \item{sicq}{The vector of Schwarz criterion with length \code{maxq} which is computed for different number of unconstrained parameters.} The following components are obtained only when \code{fixedq} is TRUE:
#' \item{muhat}{The estimator of constant parameter.}
#' \item{bhat}{The vector of evaluated constrained functional regression coefficient.}
#' \item{ahat}{The vector of unconstrained regression coefficcient estimators.}
#' \item{etahat}{The vector containing both \code{bhat} and \code{ahat} with unevaluated form.}
#' \item{yhat}{The vector of estimated response variable.}
#' \item{sp}{The vector of two penalty parameters estimated by minimising generalised cross validation (GCV).}
#' \item{L}{The number of b-spline bases used for constrained regression parameters.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{sic.l}}, \code{\link{predict.srp.c}}, \code{\link{srp.c}}
#' @examples
#' library(fda)
#' x <- matrix(rnorm(10000), ncol=100)
#' y <- matrix(rnorm(100), ncol=1)
#' p <- dim(x)[1] + 1
#' t <- seq(0, 1, length.out=dim(x)[1])*(dim(x)[1])
#' x.basis <- as.fd(splinefun(t, x[, 1], method="natural"))$basis
#' B.basis <- create.bspline.basis(rangeval=c(0, dim(x)[1]), norder=4, nbasis=35)
#' result <- sic.c(x.basis=x.basis, B.basis=B.basis, x=x, y=y, cf0=ncs(x)$cf0, maxq=10, L=35)
#' plot(result$sicq, type="b")
#' @export

sic.c <- function(x.basis=x.basis, B.basis=B.basis, x=x, y=y, cf0=cf0, maxq=maxq, fixedq=F, L=L){

  cxcf <- cf0 - rowMeans(cf0) # center the coefficients matrix
  cy <- y - mean(y)
  cx <- x - apply(x, 1, mean)
  p <- dim(x)[1]+1
  n <- dim(x)[2]

  if(fixedq==F){
    sicq <- rep(NA, maxq)

    for(i in 1:maxq){
      Vsize <- c(1:maxq)[i]
      t1 <- (seq(0, 1, length.out=p-1)*(p-1))[c(1:(p-1-Vsize))]
      J0 <- fda::inprod(x.basis, B.basis, rng=c(0, (p-1)-Vsize))
      R <- fda::eval.penalty(B.basis, fda::int2Lfd(2))
      W1 <- cbind(t(cx[c((p-Vsize):(p-1)), , drop=F]), (t(cxcf) %*% J0))
      R.s3 <- diag(Vsize + L) #scalar variables first (dim = Vsize x Vsize)

      mgfit <- mgcv::magic(y = c(cy), X = W1, sp = c(1, 1), S = list(R.s3, R), off=c(1, Vsize+1), gcv=T)

      ### Obtain the estimator of regression coefficients
      etahat <- mgfit$b
      muhat <- mean(y)-c(matrix(rowMeans(cf0), nrow=1) %*% J0 %*% matrix((etahat[-c(1:Vsize)]), ncol=1)) - c(matrix(etahat[1:Vsize], nrow=1) %*% matrix(apply(x[c((p-Vsize):(p-1)),,drop=F], 1, mean), ncol=1))

      ### SIC criterion
      mse <- mean((y-((t(cf0)%*%J0) %*% (etahat[-c(1:Vsize)])+(t(x[c((p-Vsize):(p-1)),,drop=F]) %*% matrix(etahat[1:Vsize], ncol=1)) + muhat))^2)
      sicq[i] <- log(mse)*n + log(n)*(Vsize+L+1)
    }

    qhat <- which.min(sicq)

    results <- list("qhat" = qhat, "sicq" = sicq)
    return(results)
  } else {
    qhat <- maxq
    t1 <- (seq(0, 1, length.out=p-1)*(p-1))[c(1:((p-1)-qhat))]
    J0 <- fda::inprod(x.basis, B.basis, rng=c(0, (p-1)-qhat))
    R <- fda::eval.penalty(B.basis, fda::int2Lfd(2))
    W1 <- cbind(t(cx[c((p-qhat):(p-1)), , drop=F]), (t(cf0 - rowMeans(cf0)) %*% J0))
    R.s3 <- diag(qhat + L) #scalar variables first (dim = qhat x qhat)
    mgfit <- mgcv::magic(y = c(y - mean(y)), X = W1, sp = c(1, 1), S = list(R.s3, R), off=c(1, qhat+1), gcv=T)

    ### Obtain the estimator of regression coefficients
    etahat <- mgfit$b
    bhat <- fda::eval.basis(t1, B.basis) %*% etahat[-c(1:qhat)]
    ahat <- etahat[1:qhat]
    muhat <- mean(y)-c(matrix(rowMeans(cf0), nrow=1) %*% J0 %*% matrix((etahat[-c(1:qhat)]), ncol=1)) - c(matrix(ahat, nrow=1) %*% matrix(apply(x[c((p-qhat):(p-1)),,drop=F], 1, mean), ncol=1))
    yhat <- W1%*%etahat + muhat

    results <- list("muhat" = muhat, "bhat" = bhat, "ahat" = ahat, "etahat" = etahat, "yhat" = yhat, "sp" = mgfit$sp, "L" = L)
    return(results)

  }

}
