#' Mean curve by pixel
#'
#' @param id vector of ids
#' @param year vector of year index
#' @param doy vector of day of year index
#' @param y vector of response value
#' @param rangeval a vector of length 2 containing the initial and final values of 
#' the interval over which the functional data object can be evaluated.
#' @param nbasis positive odd integer, see `create.fourier.basis` for detail
#' @param doy.eval vector of day of year index on which the mean curve is estimated. Default is 1:365
#' @param plot whether to show the plot with all trajections and mean curve, default is TRUE
#' @return a list of two elements, one is `df.mean`, a dataframe with columns `doy` and
#' `ymean`; the other is `ids.outlier`, a vector of outlier image ids
#' @export
#'
#' @examples meanCurve(tmpdf$id, tmpdf$year, tmpdf$doy, tmpdf$f1)
meanCurve <- function(x, y, x.eval, plot = FALSE) {
  ## y = f(x)
  nonna.idx = !is.na(y)
  splfit <- smooth.spline(x[nonna.idx], y[nonna.idx])
  predict(splfit, x.eval)$y
}
