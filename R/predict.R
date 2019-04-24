#' The Smooth-Rough Partition model prediction
#'
#' This function performs the predictions from the results of Smooth-Rough Partition fitting.
#'
#' The Smooth-Rough Partition model is described in "Regularised forecasting via smooth-rough partitioning of the regression coefficients", H. Maeng and P. Fryzlewicz (2018), preprint.
#'
#' @param object An object of class either 'srp.c', returned by \code{\link{srp.c}}.
#' @param x A new matrix you wish to fit Smooth-Rough Partition model. The dimension of row is the number of covariates.
#' @param  ... Further parameters that can be passed to \code{\link{predict.srp.c}}.
#' @export
#' @rdname predict.srp.c
#' @return
#' \item{yhat}{The vector of predicted values.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{sic.c}}, \code{\link{srp.c}}

predict.srp.c <- function(object, x, ...){
  Vsize <- object$qhat
  p <- dim(x)[1] + 1
  cf0 <- ncs(x)$cf0

  x.basis <- fda::as.fd(stats::splinefun((seq(0, 1, length.out=p-1)*(p-1)), x[, 1], method="natural"))$basis
  B.basis <- fda::create.bspline.basis(rangeval=c(0, p-1), norder=object$norder, nbasis=object$L)
  J0 <- fda::inprod(x.basis, B.basis, rng=c(0, (p-1)-Vsize))

  W2 <- cbind(t(x[(p-Vsize):(p-1),,drop=F]), (t(cf0)%*%J0))
  yhat <- W2%*%c(object$etahat) + object$muhat

  return(list("yhat"=yhat))
}

#' The Smooth-Rough Partition model prediction
#'
#' This function performs the predictions from the results of Smooth-Rough Partition fitting.
#'
#' The Smooth-Rough Partition model is described in "Regularised forecasting via smooth-rough partitioning of the regression coefficients", H. Maeng and P. Fryzlewicz (2018), preprint.
#'
#' @param object An object of class 'srp.l', returned by \code{\link{srp.l}}.
#' @param x A new matrix you wish to fit Smooth-Rough Partition model. The dimension of row is the number of covariates.
#' @param  ... Further parameters that can be passed to \code{\link{predict.srp.l}}.
#' @export
#' @rdname predict.srp.l
#' @return
#' \item{yhat}{The vector of predicted values.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{sic.l}}, \code{\link{srp.l}}

predict.srp.l <- function(object, x, ...){
  Vsize <- object$qhat
  p <- dim(x)[1] + 1
  cf0 <- ncs(x)$cf0

  x.basis <- fda::as.fd(stats::splinefun((seq(0, 1, length.out=p-1)*(p-1)), x[, 1], method="natural"))$basis
  m.basis = fda::create.monomial.basis(rangeval=c(0, p-1), nbasis=2)
  J0 <- fda::inprod(x.basis, m.basis, rng=c(0, (p-1)-Vsize))

  W2 <- cbind(t(x[(p-Vsize):(p-1),,drop=F]), (t(cf0)%*%J0))
  yhat <- W2%*%c(object$etahat) + object$muhat

  return(list("yhat"=yhat))
}
