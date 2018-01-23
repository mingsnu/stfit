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
#' @return a list
#' 
#' @export
#' @examples 
gapfill <- function(year, doy, mat, img.nrow, img.ncol, h,
                    doyrange = seq(min(doy), max(doy)), nnr,
                    pve = 0.99, mc.cores = parallel::detectCores()){
  idx = 1:length(year) ## idx is the index of image, 1, 2, 3,...
  registerDoParallel(cores = mc.cores)
  N = ncol(mat) # number of pixels for each image

  ## remove images that have 100% missing pixels
  pct_missing = apply(mat, 1, function(x) {sum(is.na(x))/N})
  ## Index
  idx1 = idx[pct_missing == 1] ## all missing images indexes;
  idx1c = setdiff(idx, idx1) ## not all missing images indexes;
  idx0 = idx1c[outlier(mat[idx1c,])] ## outlier images indexes;
  ## partially missing images indexes after removing outliers;
  idx2 = setdiff(idx[(pct_missing > 0) & (pct_missing < 1)], idx0) 
  ## no missing images indexes after removing outliers;
  idx3 = setdiff(idx[pct_missing == 0], idx0)
  idx4 = c(idx2, idx3)
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

  if (length(idx0) > 0)
    message(paste0(
      length(idx0),
      " outlier images are removed in the model fitting procedure."
    ))
  
  ###########################################################
  ######### Pixel-wise Temporal trend estimation ############
  ###########################################################
  ## Estimate the mean curves for each pixel
  cat("Estimating mean curve for each pixel...\n")
  ## using fully observed + partially observed - outlier images for mean estimation
  mean.mat = foreach(i = 1:N) %dopar% {
    meanCurve(doy[idx4], mat[idx4, i], doyrange)
  }
  ## mean.mat: columns are pixel index, rows are doy index (ex. 365 x 961)
  mean.mat = do.call("cbind", mean.mat)

  resid.mat = mat[idx4,]
  for(i in 1:nrow(resid.mat)){
    resid.mat[i,] = resid.mat[i,] - mean.mat[which(doyrange == doy[idx4][i]),]
  }
  
  cat("Gapfill for partial missing images...\n")
  
  ## Calculate the weight matrix
  wmat = weightMatrix(h)

  ###################################################
  ######### Covariance matrix estimation ############
  ###################################################
  cat("Estimating the covariance matrix...\n")
  ## using fully observed + partially observed - outlier images for covariance estimation
  covest = sparse_cov_est(resid.mat, img.nrow, img.ncol, nnr)
  scovest = sparseMatrix(covest$ridx, covest$cidx, x = covest$value, 
                         dims = c(N, N), symmetric = TRUE)
  varest = apply(resid.mat, 2, var, na.rm=TRUE)
  sigma2 = max(0, mean(varest - Matrix::diag(scovest))) # This is the estimated sigma_u^2 for the white noise associated with W
  
  ## Eigen decomposition of the covariance matrix; may be improved later using rARPACK
  ev = eigen(scovest)
  ev$values[ev$values < 0] = 0
  
  ev.idx = which.min(cumsum(ev$values)/sum(ev$values) < pve)
  ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
  ev.val = ev$values[1:ev.idx]
  
  ###############################################
  ######### Missing value imputation ############
  ###############################################
  ## The following uses the sparse pca to impute the missing values.
  cat("Estimating the principal component scores...\n")
  
  partial_imputed = PACE(resid.mat[1:length(idx2),], ev.vec, sigma2, ev.val, mc.cores)
  
  for(i in 1:length(idx2)){
    miss.idx = is.na(mat[idx2[i],])
    mat[idx2[i], miss.idx] = partial_imputed[i, miss.idx] + 
      mean.mat[which(doyrange == doy[idx2][i]), miss.idx]
  }
  
  # imputed_all_missing = data.frame(pixel = 1:N, df1_mean$ymean + resid_mean)
  list(imputed.partial = mat[idx2,],
       ids = list(ids.allmissing = idx1,
                  ids.partialmissing = idx2,
                  ids.fullyobserved = idx3,
                  ids.outlier = idx0),
       temporal.mean = mean.mat)
}






