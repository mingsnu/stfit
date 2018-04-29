## Perform the gapfilling algorithm
#' Gapfill missing values for a sequence of images
#'
#' @param doy vector of day of year (DOY) index
#' @param mat data matrix. Each row is a row stacked image
#' @param doyeval vector of doy to evaluate on
#' @param outlier.tol the tolerance value in defining an image as outlier. The percent of 
#' outlier pixels in an image exceed this value is regarded as outlier image which will not
#' be used in temporal mean estimation.
#' @param msk optional logistic matrix. TRUE for mask values.
#' @param cluster optional matrix defining clusters of pixels. If NULL, temporal mean estimation
#' is conducted on each pixel, otherwise all pixels from the same cluster are combined for
#' temporal mean estimation.
#' @param minimum.num.obs minimum number of observations needed for doing mean estimation for each pixel
#' @param redo whether to recalculate the mean estimation if there is an outlier.
#' @param clipRange vector of length 2, specifying the minimum and maximum values of the prediction value
#' @param clipMethod "nnr" or "truncate". "nnr" uses average of nearest neighbor pixels to impute;
#' "truncate use the clipRange value to truncate.
#' @param img.nrow number of rows for an image, only used when 'clipMethod' is "nnr"
#' @param img.ncol number of columns for an image, only used when 'clipMethod' is "nnr"
#' 
#' @return a list containing the following components:
#' `doyeval` same as input `doyeval`
#' `meanmat` matrix with number of rows equal length of `doyeval` and number of columns equal `ncol(mat)`
#' 
#' @export
#' @examples 
meanEst <- function(doy, mat,
                    doyeval = seq(min(doy), max(doy)), 
                    msk = rep(FALSE, ncol(mat)), outlier.tol = 0.5, minimum.num.obs = 4,
                    cluster = NULL, redo = TRUE, clipRange = c(-Inf, Inf), clipMethod = c("truncate", "nnr"),
                    img.nrow=NULL, img.ncol=NULL){
  idx = 1:length(doy) ## idx is the index of image, 1, 2, 3,...
  temporal_mean_est = Gapfill::opts$get("temporal_mean_est")
  N = ncol(mat) ## number of pixels
  M = nrow(mat) ## number of images
  if(length(doy) != M)
    stop("doy and mat dimension do not match!")
  if(length(msk) != N)
    stop("msk and mat dimension do not match!")
  ################################################################
  ######### Cluster/Pixel-wise overall mean estimation ###########
  ################################################################
  if(is.null(cluster)){
    ## Estimate the overall mean curves for each pixel
    cat("Estimating the overall mean curve for each pixel...\n")
    mean.mat = foreach(i = 1:N) %dopar% {
      if(msk[i])
        return(rep(NA, M)) else
          return(temporal_mean_est(doy, mat[, i], doyeval, minimum.num.obs))
    }
    ## mean.mat: columns are pixel index, rows are doy index (ex. 365 x 961)
    mean.mat = do.call("cbind", mean.mat)
    ## there will be a problem if some pixel has very few observation
    ## so check whether the missing pattern matches with 'msk'
    if(!all(msk == getMask(mean.mat)))
      warning(paste0("Some pixels which are not covered by mask have less than ", 
                     minimum.num.obs, " observations!"))
  } else {
    if(!is.vector(cluster))
      stop("cluster should be a vector.")
    if(length(cluster) != N)
      stop("cluster dimension is not correct.")
    uc = unique(cluster[!msk]) ## unique clusters w/o msk cluster
    ## Estimate the mean curves for each cluster
    cat("Estimating mean curves for each cluster...\n")
    mean.mat = matrix(NA, length(doyeval), N)
    for(cl in uc){
      clidx = cluster == cl
      mean.mat[, clidx] = temporal_mean_est(rep(doy, sum(clidx)), 
                                            c(mat[, clidx]), doyeval)
    }
  }
  
  ## residual matrix
  resid.mat = mat - mean.mat[unlist(lapply(doy, function(x,y) which(y == x), y = doyeval)),]

  ###########################################################
  ######### Outlier detection using the residuals ###########
  ###########################################################
  outlier.res = outlier(resid.mat)
  oidx = outlier.res$outpct > outlier.tol
  ## outlier image indexes
  idx0 = outlier.res$outidx[oidx]
  # outimgpct = outlier.res$outpct[oidx]
  # ooutlst = outlier.res$outlst[oidx]

  pct_missing = apply(mat, 1, function(x) {sum(is.na(x))/N})
  ## Index
  idx1 = idx[pct_missing == 1] ## all missing images indexes;
  ## partially missing image indexes after removing outliers;
  idx2 = setdiff(idx[(pct_missing > 0) & (pct_missing < 1)], idx0)
  ## no missing images indexes after removing outliers;
  idx3 = setdiff(idx[pct_missing == 0], idx0)
  
  ## Outliers are relaced with NA
  cat("Outlier pixels are replaced with NAs...\n")
  ## use NA instead of the outlier pixel values both in original data
  for(i in 1:length(outlier.res$outidx)){
    mat[outlier.res$outidx[i], outlier.res$outlst[[i]]] = NA
  }
  
  ##########################################
  ######### Redo mean estimation ###########
  ##########################################
  ## Redo pixel-wise temporal trend estimation after removing outliers
  if(redo){
    if(length(outlier.res$outidx) > 0){
      if(is.null(cluster)){
        cat("Re-estimating mean curve for each pixel...\n")
        ## using fully observed + partially observed - outlier images for mean estimation
        mean.mat = foreach(i = 1:N) %dopar% {
          if(msk[i])
            return(rep(NA, M)) else{
              if(length(idx0) > 0)
                return(temporal_mean_est(doy[-idx0], mat[-idx0, i], doyeval, minimum.num.obs)) else
                  return(temporal_mean_est(doy, mat[, i], doyeval, minimum.num.obs))
            }
        }
        mean.mat = do.call("cbind", mean.mat)
      }else{
        ## Estimate the mean curves for each pixel
        cat("Re-estimating mean curves for each cluster...\n")
        ## using fully observed + partially observed images for mean estimation
        for(cl in uc){
          clidx = cluster == cl
          if(length(idx0) > 0)
            mean.mat[, clidx] = temporal_mean_est(rep(doy[-idx0], sum(clidx)), c(mat[-idx0, clidx]), doyeval) else
              mean.mat[, clidx] = temporal_mean_est(rep(doy, sum(clidx)), c(mat[, clidx]), doyeval)
        }
      }
    }
  }
  
  ###### Correct values that beyonds clipRange
  clipMethod = match.arg(clipMethod)
  if(clipMethod == "truncate"){
    mean.mat[mean.mat < clipRange[1]] = clipRange[1]
    mean.mat[mean.mat > clipRange[2]] = clipRange[2]
  } else
    if(clipMethod == "nnr"){
      mean.mat[mean.mat < clipRange[1] | mean.mat > clipRange[2]] = NA
      ## find columns that have NA values
      msk.idx = which(msk == 1)
      notmsk.idx = which(msk == 0)
      napixel.idx = which(apply(mean.mat, 2, function(x) any(is.na(x))))
      napixel.in.msk = napixel.idx %in% msk.idx
      napixel.idx = napixel.idx[!napixel.in.msk]
      if(length(napixel.idx) > 0){
        for(i in napixel.idx){
          d = 3
          ## neighbor pixel indexes
          if(is.null(img.nrow) | is.null(img.ncol))
            stop("img.nrow and img.ncol is not supplied, which are needed for clipMethod == nnr.")
          nbrpixel.idx = intersect(nbr(i-1, img.nrow, img.ncol, d, d) + 1, notmsk.idx)
          while(all(nbrpixel.idx %in% napixel.idx)){
            d = d + 2
            nbrpixel.idx = intersect(nbr(i-1, img.nrow, img.ncol, d, d) + 1, notmsk.idx)
          }
          ## mean of neighborhood pixels
          mm = apply(mean.mat[, nbrpixel.idx], 1, mean, na.rm=TRUE)
          mm.miss.idx = is.na(mean.mat[,i])
          mean.mat[mm.miss.idx, i] = mm[mm.miss.idx]
        }
      }
    }
  

  if (length(idx0) > 0)
    message(paste0(
      length(idx0),
      " outlier images are removed in the mean estimation procedure."))
  return(list(
    doyeval = doyeval,
    meanmat = mean.mat,
    idx = list(idx.allmissing = idx1,
               idx.partialmissing = idx2,
               idx.fullyobserved = idx3,
               idx.outlier = idx0),
    outlier = outlier.res
  ))
}




