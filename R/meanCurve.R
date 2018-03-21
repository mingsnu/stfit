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
  ## TODO: x.eval should be witin certain range.
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > 4){
    splfit <- smooth.spline(x[nonna.idx], y[nonna.idx])
    res = predict(splfit, x.eval)$y
    # res[x.eval < min(x) | x.eval > max(x)] = NA
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}