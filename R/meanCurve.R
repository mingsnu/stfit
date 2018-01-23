#' Mean curve by pixel
#'
#' @param x independent variable
#' @param y dependent variable
#' @param x.eval vector to predict on
#' @param plot logical. NOT IMPLEMENTED YET
#'
#' @export
#'
#' @examples 
meanCurve <- function(x, y, x.eval, plot = FALSE) {
  ## y = f(x)
  nonna.idx = !is.na(y)
  splfit <- smooth.spline(x[nonna.idx], y[nonna.idx])
  predict(splfit, x.eval)$y
}
