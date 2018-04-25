#' Smoothing spline regression
#'
#' @param x independent variable
#' @param y response variable
#' @param x.eval vector to predict on
#' @param plot logical. NOT IMPLEMENTED YET
#'
#' @return predicted values at 'x.eval'
#' @export
#'
#' @examples 
smooth_spline <- function(x, y, x.eval=x, minimum.num.obs=4, ...) {
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > minimum.num.obs){
    splfit <- smooth.spline(x[nonna.idx], y[nonna.idx], ...)
    res = predict(splfit, x.eval)$y
    # res[x.eval < min(x) | x.eval > max(x)] = NA
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

#' Local linear regression
#'
#' @param x independent variable
#' @param y response variable
#' @param h bandwidth
#' @param Kern Kernel
#' @param x.eval dnew data to predict on
#'
#' @return predicted values at 'x.eval'
#' @export
#'
#' @examples
#' 
llreg <- function(x, y, x.eval=x, minimum.num.obs = 4, h=60, Kern=epan){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > minimum.num.obs){
    res <- .llreg(x[nonna.idx], y[nonna.idx], x.eval, h, Kern)
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

.llreg <- function(x, y, x.eval, h, Kern){
  if(length(x) != length(y))
    stop("llreg: x and y should have the same length")
  N = length(x.eval)
  fitted = rep(0, N)
  for(i in 1:N){
    vecX = x - x.eval[i]
    vecK = Kern((vecX)/h)
    sumK = sum(vecK)
    h1 = h
    while(sumK == 0){
      h1 = h1 + 1
      vecK = Kern((vecX)/h1)
      sumK = sum(vecK)
    }
    S0 = mean(vecK)
    S1 = mean(vecX*vecK) ## higher order term than S0 and S2, keep this term for accuracy 
    S2 = mean(vecX^2*vecK)
    denom = S2*S0 - S1^2
    if(denom != 0){
      fitted[i] = mean((S2-S1*vecX)*vecK*y)/denom ## local linear estimator
    } else
      fitted[i] = sum(vecK*y)/sumK ## local constant estimator
  }
  fitted
}

#' Local Polynomial Regression
#'
#' @param x independent variable
#' @param y response variable
#' @param x.eval vector to predict on
#' @param span see 'loess' function
#' @param ... other parameters passed to 'loess' function
#'
#' @return predicted values at 'x.eval'
#' @export
#'
#' @examples
lpreg <- function(x, y, x.eval, minimum.num.obs = 4, span=0.3, ...){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > minimum.num.obs){
    x = x[nonna.idx]
    y = y[nonna.idx]
    loessfit <- loess(y~x, span = span, control = loess.control(surface = "direct"), ...)
    res = predict(loessfit, data.frame(x = x.eval))
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

#' spline regression
#'
#' @param x independent variable
#' @param y response variable
#' @param x.eval vector to predict on
#' @param rangeeval see \code{fda::create.basis}
#' @param nbasis see \code{fda::create.basis}
#' @param ... arguments passed to \code{fad::create.basis} functions
#'
#' @return predicted values at 'x.eval'
#' @export
#'
#' @examples
spreg <- function(x, y, x.eval, minimum.num.obs = 4, basis = c("fourier", "bspline"), 
                  rangeval = c(min(x.eval)-1, max(x.eval)), nbasis = 11, ...){
  nonna.idx = !is.na(y)
  basis = match.arg(basis)
  if(sum(nonna.idx) > minimum.num.obs){
    x = x[nonna.idx]
    y = y[nonna.idx]
    # whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
    # idx = (y > whisker[1]) & (y <whisker[2])
    # x = x[idx]
    # y = y[idx]
    if(basis == "fourier"){
      bs = fda::create.fourier.basis(rangeval=rangeval, nbasis=nbasis, ...)
    } else
      if(basis == "bspline"){
        bs = fda::create.bspline.basis(rangeval=rangeval, nbasis=nbasis, ...)
      }
    X = fda::eval.basis(x, bs)
    lmfit = lm.fit(X, y)
    return(fda::eval.basis(x.eval, bs) %*% lmfit$coefficients)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

