#' Local constant estimation for image mean
#'
#' @param mat a matrix, each column is a stacked values of an image
#' @param wmat weight matrix
#'
#' @return
#' @export
#'
#' @examples
meanEst <- function(mat, wmat) {
  mean_est(mat, wmat)
}