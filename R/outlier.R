## returns outlier image index (row index of mat)
#' Title
#'
#' @param mat data matrix. Each row is a row stacked image.
#'
#' @return row indexes for outlier images
#' @export
#'
#' @examples
outlier <- function(mat){
  ## If any pixel of an image is an outlier, the image is treated
  ## as an outlier image and is removed from further calculation
  .outlier <- function(y){
    whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
    which(y < whisker[1] | y > whisker[2])
  }
  sort(unique(unlist(lapply(1:ncol(mat), function(i) .outlier(mat[,i])))))
}
