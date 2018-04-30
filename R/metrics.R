#' Error metrics
#'
#' @param y vector
#' @param ypred vecotr
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#'
#' @examples
RMSE = function(y, ypred){
  sqrt(sum((y-ypred)^2)/length(y))
}

#' Error metrics
#'
#' @param y vector
#' @param ypred vector
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#'
#' @examples
NMSE = function(y, ypred){
  sum((y-ypred)^2)/sum(y^2)
}

#' Title
#'
#' @param y vector
#' @param ypred vector
#'
#' @return numeric number. A measure of difference between y and ypred.
#' @export
#'
#' @examples
ARE = function(y, ypred){
  mean(abs(y - ypred)/y)
}