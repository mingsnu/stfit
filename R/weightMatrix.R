#' weight matrix calculation
#'
#' @param h 'bandwith'
#'
#' @return a matrix indecating the neiborhood weighting structure
#'
#' @examples
weightMatrix <- function(h){
  matrix(c(0.375, 0.562, 0.375, 0.562, 0.75, 0.562, 0.375, 0.562, 0.375), 3)
}