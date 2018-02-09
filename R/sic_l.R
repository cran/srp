#' The smoothness change-point detection of regression coefficients in (simpler) Smooth-Rough Partition linear regression model
#'
#' This function performs the optimisation of the number of unconstrained regression parameters in (simpler) Smooth-Rough Partition model and gives the change-point in regression parameters.
#'
#' The Smooth-Rough Partition model is described in "Regularised forecasting via smooth-rough partitioning of the regression coefficients", H. Maeng and P. Fryzlewicz (2018), preprint.
#'
#' @param x.basis The b-spline basis defined for interpolated x in \code{\link{srp.l}}.
#' @param M.basis The monomial basis defined for constrained regression coefficient.
#' @param x The design matrix used in \code{\link{srp.l}}.
#' @param y The response variable used in \code{\link{srp.l}}.
#' @param cf0 The coefficient matrix obtained by natural cubic spline interpolation of x in \code{\link{ncs}}.
#' @param maxq The maximum value of possible number of unconstrained parameters if \code{fixedq} is FALSE. Otherwise, it is considered as an unique number of unconstrained parameters.
#' @param fixedq If TRUE, \code{maxq} is considered as a fixed number of unconstrained parameters and if FALSE, \code{maxq} is a maximum and a sequence of possible values are investigated to select the optimal.
#' @return The following components are obtained only when \code{fixedq} is FALSE:
#' \item{qhat}{The optimal number of unconstrained parameters.}
#' \item{sicq}{The vector of Schwarz criterion with length \code{maxq} which is computed for different number of unconstrained parameters.} The following components are obtained only when \code{fixedq} is TRUE:
#' \item{muhat}{The estimator of constant parameter.}
#' \item{bhat}{The vector of evaluated constrained functional regression coefficient.}
#' \item{ahat}{The vector of unconstrained regression coefficcient estimators.}
#' \item{etahat}{The vector containing both \code{bhat} and \code{ahat} with unevaluated form.}
#' \item{yhat}{The vector of estimated response variable.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{sic.c}}, \code{\link{srp.l}}
#' @examples
#' library(fda)
#' x <- matrix(rnorm(10000), ncol=100)
#' y <- matrix(rnorm(100), ncol=1)
#' p <- dim(x)[1] + 1
#' t <- seq(0, 1, length.out=dim(x)[1])*(dim(x)[1])
#' x.basis <- as.fd(splinefun(t, x[, 1], method="natural"))$basis
#' M.basis <- create.monomial.basis(rangeval=c(0, dim(x)[1]), nbasis=2)
#' result <- sic.l(x.basis=x.basis, M.basis=M.basis, x=x, y=y, cf0=ncs(x)$cf0, maxq=10)
#' plot(result$sicq, type="b")
#' @export

sic.l <- function(x.basis=x.basis, M.basis=M.basis, x=x, y=y, cf0=cf0, maxq=maxq, fixedq=F){

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
      J0 <- fda::inprod(x.basis, M.basis, rng=c(0, (p-1)-Vsize))
      W1 <- cbind(t(cx[c((p-Vsize):(p-1)), , drop=F]), (t(cxcf)%*%J0))

      ### Obtain the estimator of regression coefficients
      etahat <- solve(t(W1) %*% W1) %*% t(W1) %*% cy
      bhat <- fda::eval.basis(t1, M.basis)%*%etahat[-c(1:Vsize)]
      ahat <- etahat[1:Vsize]
      ahatv <- c(matrix(ahat, nrow=1) %*% matrix(apply(x[c((p-Vsize):(p-1)),,drop=F], 1, mean), ncol=1)) # t(alpha hat) %*% vbar
      bhatx <- c(matrix(rowMeans(cf0), nrow=1) %*% J0 %*% matrix(etahat[-c(1:Vsize)], ncol=1))
      muhat <- mean(y) - ahatv - bhatx

      ### SIC criterion
      sicq[i] = log(mean((y - (W1%*%etahat + muhat))^2))*n + log(n)*(Vsize+2+1)
    }

    qhat <- which.min(sicq)

    results <- list("qhat" = qhat, "sicq" = sicq)
    return(results)
  } else {
    qhat <- maxq
    t1 <- (seq(0, 1, length.out=p-1)*(p-1))[c(1:((p-1)-qhat))]
    J0 <- fda::inprod(x.basis, M.basis, rng=c(0, (p-1)-qhat))
    W1 <- cbind(t(cx[c((p-qhat):(p-1)), , drop=F]), (t(cf0 - rowMeans(cf0)) %*% J0))

    ### Obtain the estimator of regression coefficients
    etahat <- solve(t(W1) %*% W1) %*% t(W1) %*% (y - mean(y))
    bhat <- fda::eval.basis(t1, M.basis)%*%etahat[-c(1:qhat)]
    ahat <- etahat[1:qhat]
    ahatv <- c(matrix(ahat, nrow=1) %*% matrix(apply(x[c((p-qhat):(p-1)),,drop=F], 1, mean), ncol=1)) # t(alpha hat) %*% vbar
    bhatx <- c(matrix(rowMeans(cf0), nrow=1) %*% J0 %*% matrix(etahat[-c(1:qhat)], ncol=1))
    muhat <- mean(y) - ahatv - bhatx

    yhat <- W1%*%etahat + muhat # test set

    results <- list("muhat" = muhat, "bhat" = bhat, "ahat" = ahat, "etahat" = etahat, "yhat" = yhat)
    return(results)

  }

}
