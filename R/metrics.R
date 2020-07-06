#' Root Mean Square Estimation
#'
#' @param y vector
#' @param ypred vecotr
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#'
#' @examples
RMSE = function(y, ypred){
  idx = !is.na(y) & !is.na(ypred)
  sqrt(sum((y[idx]-ypred[idx])^2)/sum(idx))
}

#' Normalized Mean Square Estimation
#'
#' @param y vector
#' @param ypred vector
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#'
#' @examples
NMSE = function(y, ypred){
  idx = !is.na(y) & !is.na(ypred)
  sum((y[idx]-ypred[idx])^2)/sum(y[idx]^2)
}

#' Absolute relative error
#'
#' @param y vector
#' @param ypred vector
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#'
#' @examples
ARE = function(y, ypred){
  idx = !is.na(y) & !is.na(ypred)
  mean(abs(y[idx] - ypred[idx])/y[idx])
}