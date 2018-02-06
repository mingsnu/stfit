## returns outlier image index (row index of mat)
#' Title
#'
#' @param mat data matrix. Each row is a row stacked image.
#'
#' @return row indexes for outlier images
#' @export
#'
#' @examples
# outlier <- function(mat){
#   ## If any pixel of an image is an outlier, the image is treated
#   ## as an outlier image and is removed from further calculation
#   .outlier <- function(y){
#     whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
#     which(y < whisker[1] | y > whisker[2])
#   }
#   sort(unique(unlist(lapply(1:ncol(mat), function(i) .outlier(mat[,i])))))
# }
outlier <- function(mat, outlier.tol){
  .outlier <- function(y){
    whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
    which(y < whisker[1] | y > whisker[2])
  }
  ## a table of the outlier row index of mat, higher value indicate higher probability that the
  ## corresponding image has problem.
  lst = lapply(1:ncol(mat), function(i) .outlier(mat[,i])) ## length equals # of pixels
  tbl = table(unlist(lst))
  idx = as.numeric(names(tbl))
  tot = apply(mat[idx, ], 1, function(x) sum(!is.na(x)))
  outpct = tbl/tot
  outidx = idx[outpct > outlier.tol]
  outlst = vector("list", length(outidx))
  for(i in 1:length(outidx)){
    for(j in 1:length(lst)){
      if(outidx[i] %in% lst[[j]])
        outlst[[i]] = c(outlst[[i]], j)
    }
  }
  return(list(outidx = outidx, outpct = outpct[outpct > outlier.tol], outlst = outlst))
}
