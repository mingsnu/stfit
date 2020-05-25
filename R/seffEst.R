## Perform the gapfilling algorithm
#' Gapfill missing values for a sequence of images
#' 
#' @param rmat residual matrix
#' @param img.nrow image row dimension
#' @param img.ncol image column dimension
#' @param pve percent of variance explained of the selected eigen values
#' @param h.cov bandwidth for spatial covariance estimation; ignored if `weight.cov` is supplied
#' @param h.sigma2 bandwidth for sigma2 estimation
#' @param weight.cov weight vector for spatial covariance estimation
#' @param weight.sigma2 weight vector for sigma2 estimation;
#' @param nnr maximum of nearest neibors used to calculate correlation
#' @param method "lc" for local constant covariance estimation and "emp" for empirical covariance estimation
#' @param msk optional logistic matrix. TRUE for mask values.
#' @param msk.tol if 'msk' is not given, the program will determine the mask using 'getMask'
#' function. If the percentage of missing values for a pixel over time is greater than this
#' @param keep.original whether to keep the originally observed values in the resulting matrix
#' @param detail if TRUE returns eigen values,eigen functions
#' @return a list containing the following components:
#' `year` same as input `year`
#' `doy` same as input `doy`
#' `imputed.partial`
#' 
#' @export
seffEst <- function(rmat, img.nrow, img.ncol, h.cov = 2, h.sigma2 = 2,
                    weight.cov = NULL, weight.sigma2 = NULL,
                    nnr, method = c("lc", "emp"), keep.original = FALSE,
                    pve = 0.99, msk = NULL, msk.tol = 0.95,
                    detail = FALSE){
  cat("Estimating spatial effect...")
  if(is.null(weight.cov))
    weight.cov = weightMatrix(h.cov)
  if(is.null(weight.sigma2))
    weight.sigma2 = weightVector(h.sigma2)
  ## initialize the rmat imputed matrix
  idx = 1:nrow(rmat)
  N = img.nrow * img.ncol ## total number of pixels (including pixels in mask)
  if(ncol(rmat) != N)
    stop("number of columns of rmat does not match img.nrow*img.ncol.")

  if(is.null(msk)){
    msk = getMask(rmat, tol = msk.tol) # idx for 'black holes'
  } else {
    if(!is.vector(msk))
      stop("msk should be a vector.")
    if(length(msk) != N)
      stop("msk dimension is not correct.")
  }
  ## number of 'actual' pixels for each image (except for pixels in mask)
  N1 = sum(!msk) 
  
  ## initial spatial effect matrix
  seffmat = rmat
  seffmat[is.na(seffmat)] = 0
  
  ## no spatial effect if there are too few actual pixels in one image
  if(N1 <= 1)
    return(list(seffmat=seffmat,
                idx = list()))
  
  pidx = (1:N)[!msk] ## 'actual' pixel indexes
  ## keep 'actual' pixels only, 
  ## The i-th column in the old rmat correspones to the
  ## `which(pidx == i)`th column in the new rmat
  rmat = rmat[,!msk, drop=FALSE] ## now 'rmat' has N1 columns
  
  ## remove images that have 100% missing pixels
  pct_missing = apply(rmat, 1, function(x) {sum(is.na(x))/N1})
  ## Index
  idx1 = idx[pct_missing == 1] ## all missing images indexes;
  idx2 = idx[pct_missing > 0 & pct_missing < 1] ## partially missing images indexes;
  idx3 = idx[pct_missing == 0] ## non missing images indexes;
  
  ###################################################
  ######### Covariance matrix estimation ############
  ###################################################
  cat("Estimating the covariance matrix...\n")
  ## using fully observed + partially observed images for covariance estimation
  method <- match.arg(method)
  if(method == "lc"){
    if(N == N1)
      covest = sparse_lc_cov_est(rmat, weight.cov, img.nrow, img.ncol, nnr)
    else
      covest = sparse_lc_cov_est1(rmat, weight.cov, img.nrow, img.ncol, nnr, pidx - 1)
  } else
    if (method == "emp") {
      if (N == N1)
        ## no black hole
        covest = sparse_emp_cov_est(rmat, img.nrow, img.ncol, nnr)
      else
        covest = sparse_emp_cov_est1(rmat, img.nrow, img.ncol, nnr, pidx -
                                       1) #pidx start from 0 in C++
    }
  ## coerce NA to 0; two vectors with no overlapped observations leads to NA
  covest$value[is.na(covest$value)] = 0
  scovest = sparseMatrix(covest$ridx, covest$cidx, x = covest$value, 
                         dims = c(N1, N1), symmetric = TRUE)
  cat("Estimating the variance function...\n")
  ## TODO: using local linear smoothing
  sigma2 = apply(rmat, 2, var, na.rm=TRUE)
  nugg = max(0, mean(sigma2 - Matrix::diag(scovest), na.rm=TRUE))
  
  cat("Doing eigen decomposition on covariance matrix...\n")
  ## Eigen decomposition of the covariance matrix; may be improved later using rARPACK
  ev = eigen(scovest)
  ev$values[ev$values < 0] = 0
  
  ev.idx = which.min(cumsum(ev$values)/sum(ev$values) < pve)
  cat("Spatial effect covariance estimation: the first ", ev.idx, " eigen values are used...\n")
  ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
  ev.val = ev$values[1:ev.idx]
  
  ###############################################
  ######### Missing value imputation ############
  ###############################################
  ## The following uses the sparse pca to impute the missing values.
  cat("Estimating the principal component scores for partially missing images...\n")
  rmat_imputed_partial = PACE(rmat[idx2,, drop = FALSE], ev.vec, nugg, ev.val)
  if(keep.original){
    partial_miss_idx = is.na(rmat[idx2,,drop=FALSE])
    seffmat[idx2, !msk][partial_miss_idx] = rmat_imputed_partial[partial_miss_idx]
  } else {
    seffmat[idx2, !msk] = rmat_imputed_partial
  }

  if(detail){
      return(list(
        seffmat = seffmat,
        idx = list(idx.allmissing = idx1,
                   idx.partialmissing = idx2,
                   idx.fullyobserved = idx3),
        eigen.vec = ev.vec, 
        nugg = nugg, 
        eigen.val = ev.val))
  } else {
      return(list(
        seffmat = seffmat,
        idx = list(idx.allmissing = idx1,
                   idx.partialmissing = idx2,
                   idx.fullyobserved = idx3)
      ))
  }
}




