## Perform the gapfilling algorithm
#' Gapfill missing values for a sequence of images
#'
#' @param year vector of year index
#' @param doy vector of day of year (DOY) index
#' @param mat data matrix. Each row is a row stacked image
#' @param img.nrow image row dimension
#' @param img.ncol image column dimension
#' @param h bandwidth
#' @param doyrange vector of two elements, minimum and maximum value of DOY
#' @param pve percent of variance explained of the selected eigen values
#' @param mc.cores number of cores to use for computation
#' 
#' @return a list containing the following components:
#' `year` same as input `year`
#' `doy` same as input `doy`
#' `imputed.partial`
#' 
#' @export
#' @examples 
gapfill <- function(year, doy, mat, img.nrow, img.ncol, h,
                    doyrange = seq(min(doy), max(doy)), nnr, method = c("lc", "emp"),
                    pve = 0.99, outlier.tol = 0.2, outlier.action = c("ac", "keep"), mc.cores = parallel::detectCores()){
  idx = 1:length(year) ## idx is the index of image, 1, 2, 3,...
  registerDoParallel(cores = mc.cores)
  msk = getMask(mat) # idx for 'black holes'
  N = img.nrow * img.ncol ## total number of pixels (including black holes)
  N1 = sum(!msk) # number of 'actual' pixels for each image (except for black holes)
  if(N1 == 0){ ## a black hole image
    return(list(
      year = year, 
      doy = doy,
      imputed.mat = mat,
      idx = list(idx.allmissing = idx,
                 idx.partialmissing = c(),
                 idx.fullyobserved = c(),
                 idx.outlier = c()),
      temporal.mean = NULL))
  }
  
  pidx = (1:N)[!msk] ## 'actual' pixel indexes
  ## keep 'actual' pixels only, 
  ## The i-th column in the old mat correspones to the
  ## `which(pidx == i)`th column in the new mat
  mat = mat[,!msk, drop=FALSE] ## now 'mat' has N1 columns

  ## remove images that have 100% missing pixels
  pct_missing = apply(mat, 1, function(x) {sum(is.na(x))/N1})
  ## Index
  idx1 = idx[pct_missing == 1] ## all missing images indexes;
  idx1c = setdiff(idx, idx1) ## not all missing images indexes;
  if (length(idx1) > 0)
    message(
      paste0(
        length(idx1),
        " out of ",
        length(pct_missing),
        " images have no observations and are removed in the model fitting procedure;",
        " global means are used in imputing these images."
      )
    )
  
  ###########################################################
  ######### Pixel-wise Temporal trend estimation ############
  ###########################################################
  ## Estimate the mean curves for each pixel
  cat("Estimating mean curve for each pixel...\n")
  ## using fully observed + partially observed images for mean estimation
  mean.mat = foreach(i = 1:N1) %dopar% {
    meanCurve(doy[idx1c], mat[idx1c, i], doyrange)
  }
  ## mean.mat: columns are pixel index, rows are doy index (ex. 365 x 961)
  mean.mat = do.call("cbind", mean.mat)
  ## find columns that have NA values
  napixel.idx = which(apply(mean.mat, 2, function(x) any(is.na(x))))
  while(length(napixel.idx) > 0){
    for(i in napixel.idx){
      d = 3
      ## neighbor pixel indexes
      nbrpixel.idx = intersect(nbr(pidx[i]-1, img.nrow, img.ncol, d, d) + 1, pidx)
      col.idx = which(pidx %in% nbrpixel.idx)
      ## mean of neighborhood pixels
      mm = apply(mean.mat[,col.idx], 1, FUN = function(x){
        if(all(is.na(x)))
          return(NA) else
            return(mean(x, na.rm = TRUE))
      })
      mm.miss.idx = is.na(mean.mat[,i])
      mean.mat[mm.miss.idx, i] = mm[mm.miss.idx]
    }
    napixel.idx = which(apply(mean.mat, 2, function(x) any(is.na(x))))
  }
  ## residual matrix
  resid.mat = mat[idx1c,]
  for(i in 1:nrow(resid.mat)){
    resid.mat[i,] = resid.mat[i,] - mean.mat[which(doyrange == doy[idx1c][i]),]
  }
  ###########################################################
  ######### Outlier detection using the residuals ###########
  ###########################################################
  outlier.res = outlier(resid.mat, outlier.tol)
  ## outlier images indexes;
  if(length(outlier.res$outidx) > 0)
    idx0 = idx1c[outlier.res$outidx] else
      idx0 = c()
  ## partially missing image indexes after removing outliers;
  idx2 = setdiff(idx[(pct_missing > 0) & (pct_missing < 1)], idx0)
  ## no missing images indexes after removing outliers;
  idx3 = setdiff(idx[pct_missing == 0], idx0)
  idx4 = c(idx2, idx3)
  
  ## Redo pixel-wise temporal trend estimation after removing outliers
  if(length(outlier.res$outidx) > 0){
    cat("Re-estimating mean curve for each pixel...\n")
    ## using fully observed + partially observed - outlier images for mean estimation
    mean.mat = foreach(i = 1:N1) %dopar% {
      meanCurve(doy[idx4], mat[idx4, i], doyrange)
    }
    mean.mat = do.call("cbind", mean.mat)
    ## find columns that have NA values
    napixel.idx = which(apply(mean.mat, 2, function(x) any(is.na(x))))
    while(length(napixel.idx) > 0){
      for(i in napixel.idx){
        d = 3
        ## neighbor pixel indexes
        nbrpixel.idx = intersect(nbr(pidx[i]-1, img.nrow, img.ncol, d, d) + 1, pidx)
        col.idx = which(pidx %in% nbrpixel.idx)
        ## mean of neighborhood pixels
        mm = apply(mean.mat[,col.idx], 1, FUN = function(x){
          if(all(is.na(x)))
            return(NA) else
              return(mean(x, na.rm = TRUE))
        })
        mm.miss.idx = is.na(mean.mat[,i])
        mean.mat[mm.miss.idx, i] = mm[mm.miss.idx]
      }
      napixel.idx = which(apply(mean.mat, 2, function(x) any(is.na(x))))
    }
    resid.mat = mat[idx4,]
    for(i in 1:nrow(resid.mat)){
      resid.mat[i,] = resid.mat[i,] - mean.mat[which(doyrange == doy[idx4][i]),]
    }
  }

  if (length(idx0) > 0)
    message(paste0(
      length(idx0),
      " outlier images are removed in the model fitting procedure."
    ))
  
  cat("Gapfill for partial missing images...\n")
  
  ## Calculate the weight matrix
  wmat = weightMatrix(h)
  
  ###################################################
  ######### Covariance matrix estimation ############
  ###################################################
  cat("Estimating the covariance matrix...\n")
  ## using fully observed + partially observed - outlier images for covariance estimation
  method <- match.arg(method)
  if(method == "emp"){
    if(N == N1) ## no black hole
      covest = sparse_emp_cov_est(resid.mat, img.nrow, img.ncol, nnr) else
        covest = sparse_emp_cov_est1(resid.mat, img.nrow, img.ncol, nnr, pidx-1) #pidx start from 0 in C++
  }else
      if(method == "lc")
        if(N == N1)
          covest = sparse_lc_cov_est(resid.mat, wmat, img.nrow, img.ncol, nnr) else
            covest = sparse_lc_cov_est1(resid.mat, wmat, img.nrow, img.ncol, nnr, pidx - 1)
        
  scovest = sparseMatrix(covest$ridx, covest$cidx, x = covest$value, 
                         dims = c(N1, N1), symmetric = TRUE)
  cat("Estimating the variance function...\n")
  varest = apply(resid.mat, 2, var, na.rm=TRUE)
  sigma2 = max(0, mean(varest - Matrix::diag(scovest))) # This is the estimated sigma_u^2 for the white noise associated with W
  
  cat("Doing eigen decomposition on covariance matrix...\n")
  ## Eigen decomposition of the covariance matrix; may be improved later using rARPACK
  ev = eigen(scovest)
  ev$values[ev$values < 0] = 0
  
  ev.idx = which.min(cumsum(ev$values)/sum(ev$values) < pve)
  cat("The first ", ev.idx, " eigen values are used...\n")
  ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
  ev.val = ev$values[1:ev.idx]
  
  ###############################################
  ######### Missing value imputation ############
  ###############################################
  ## The following uses the sparse pca to impute the missing values.
  cat("Estimating the principal component scores for partially missing images...\n")
  partial_imputed = PACE(resid.mat[1:length(idx2),, drop = FALSE], ev.vec, sigma2, ev.val)
  cat("Gapfilling partially missing images...\n")
  for(i in 1:length(idx2)){
    miss.idx = is.na(mat[idx2[i],])
    mat[idx2[i], miss.idx] = partial_imputed[i, miss.idx] + 
      mean.mat[which(doyrange == doy[idx2][i]), miss.idx]
  }
  if(length(outlier.res$outidx) > 0){
    outlier.action <- match.arg(outlier.action)
    outlier.resid.mat = mat[idx0,, drop=FALSE]
    if(outlier.action == "ac"){
      cat("Outlier autocorrection is used, outlier pixels are treated as missing...\n")
      ## use NA instead of the outlier pixel values both in original data and residual data
      for(i in 1:length(outlier.res$outlst)){
        mat[idx0[i], outlier.res$outlst[[i]]] = NA
        outlier.resid.mat[i, outlier.res$outlst[[i]]] = NA
      }
    } else if(outlier.action == "keep"){
      cat("Outlier pixels are kept...\n")
      for(i in 1:length(outlier.res$outlst)){
        outlier.resid.mat[i, outlier.res$outlst[[i]]] = NA
      }
    } else
      stop("not defined method for outlier.action.")
    for (i in 1:nrow(outlier.resid.mat)) {
      outlier.resid.mat[i, ] = outlier.resid.mat[i, ] - mean.mat[which(doyrange == doy[idx0][i]), ]
    }
    ## if all pixels of an image are outliers, the image is imputed with mean
    tmpidx = outlier.res$outpct == 1
    if (sum(tmpidx) > 0) {
      cat("All pixels in image with doy = ", doy[idx0][tmpidx],
          " are outliers. Temporal mean surface is used to impute it.")
      # outlier.resid.mat[tmpidx,] = mean.mat[which(doyrange == doy[idx0][tmpidx]),]
      if(sum(tmpidx) != length(tmpidx)){
        outlier_imputed = matrix(0, nrow(outlier.resid.mat), ncol(outlier.resid.mat))
        cat("Estimating the principal component scores for outlier missing images...\n")
        outlier_imputed[!tmpidx, ] = PACE(outlier.resid.mat[!tmpidx,,drop=FALSE], ev.vec, sigma2, ev.val)
      } else
        outlier_imputed = matrix(0, nrow(outlier.resid.mat), ncol(outlier.resid.mat))
    } else{
      cat("Estimating the principal component scores for outlier missing images...\n")
      outlier_imputed = PACE(outlier.resid.mat, ev.vec, sigma2, ev.val)
    }
    cat("Gapfilling outlier missing images...\n")
    for(i in 1:length(idx0)){
      miss.idx = is.na(mat[idx0[i],])
      mat[idx0[i], miss.idx] = outlier_imputed[i, miss.idx] + 
        mean.mat[which(doyrange == doy[idx0][i]), miss.idx]
    }
  }
  
  if(N == N1){
    return(list(
      year = year, 
      doy = doy,
      imputed.mat = mat,
         idx = list(idx.allmissing = idx1,
                    idx.partialmissing = idx2,
                    idx.fullyobserved = idx3,
                    idx.outlier = idx0),
      outlier.lst = lapply(outlier.res$outlst, function(x) pidx[x]),
      temporal.mean = mean.mat,
      eigen.vec = ev.vec, 
      sigma2 = sigma2, 
      eigen.val = ev.val))
  } else{
    impmat = matrix(NA, nrow(mat), length(msk))
    impmat[, !msk] = mat
    
    mmat = matrix(NA, nrow(mean.mat), length(msk))
    mmat[,!msk] = mean.mat
    # imputed_all_missing = data.frame(pixel = 1:N, df1_mean$ymean + resid_mean)
    return(list(
      year = year, 
      doy = doy,
      imputed.mat = impmat,
         idx = list(idx.allmissing = idx1,
                    idx.partialmissing = idx2,
                    idx.fullyobserved = idx3,
                    idx.outlier = idx0),
      outlier.lst = lapply(outlier.res$outlst, function(x) pidx[x]),
      temporal.mean = mmat,
      eigen.vec = ev.vec, 
      sigma2 = sigma2, 
      eigen.val = ev.val))
  }
}




