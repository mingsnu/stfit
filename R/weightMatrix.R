#' weight matrix calculation
#'
#' @param h 'bandwith'
#'
#' @return a weighting matrix
weightMatrix <- function(h){
  ## matrix(c(0.375, 0.562, 0.375, 0.562, 0.75, 0.562, 0.375, 0.562, 0.375), 3)
  if(h <=0)
    stop("bandwidth can not be non-positive.")
  if(h == 1)
    warning("No neighborhood information is used")
  aa = seq(-h, h, by = 1)/h
  out = outer(aa, aa, function(x, y) epan(sqrt(x^2+y^2)))
  return(out[2:(2*h), 2:(2*h)])
}

#' weight vector calculation
#'
#' @param h bandwidth, should be positive numbers
#'
#' @return a vector
#' @export
weightVector <- function(h){
  if(h <=0)
    stop("bandwidth can not be non-positive.")
  if(h == 1)
    warning("No neighborhood information is used")
  out = epan(seq(-h, h, by = 1)/h)
  return(out[2:(2*h)])
}


