#' Gapfilling for repeated measurements over years
#'
#' @param year vector of year
#' @param doy vector of doy
#' @param mat matrix with number of row equals length of doy 
#' @param img.nrow row number for the image to gapfill
#' @param img.ncol column number for the image to gapfill
#' @param doyeval a vector of DOY on which to perform the gapfill
#' @param h.tcov bandwidth for temporal covariance estimation
#' @param h.tsigma2 bandwith for temporal variance estimation 
#' @param h.scov bandwidth for spatial covariance estimation
#' @param h.ssigma2 bandwidth for spatial variance estimation
#' @param nnr number of nearest neighbor pixels to use for spatial covariance estimation
#' @param outlier.action "keep" to keep outliers; "remove" to replace outliers with gapfilled values
#' @param intermediate.save TRUE or FASLE; whether to save the intermediate results
#' @param intermediate.dir directory where to save the intermediate results
#' @param use.intermediate.result whether to use the intermediate results in the 'intermediate.dir' folder instead
#' of recalcuating 
#' @param teff TRUE or FALSE, wheter to calculate the temporal effect
#' @param seff TRUE or FALSE, wheter to calculate the spatial effect
#' @param doy.break a vector of break points for `doy`. spatial effect are estimated seperately on each interval
#' @param b2e.con TRUE or FALSE. whether the the begining of `doy` is connected with the end of `doy`.
#'
#' @return List of length 3 with entries:
#' \itemize{
#'   \item imat imputed matrix of `mat`
#'   \item idx a list of image indexes information
#'   \item outlier a list of image outliers information
#' }
#' @export
#'
#' @examples
gapfill_landsat <- function(year, doy, mat, img.nrow, img.ncol, doyeval = 1:365,  h.tcov = 100, h.tsigma2 = 300,
                            h.scov = 2, h.ssigma2 = 2, nnr = 10, outlier.action = c("keep", "remove"), outlier.tol = 0.2,
                            intermediate.save = TRUE, intermediate.dir = "./output/",
                            use.intermediate.result = TRUE, teff = TRUE, seff = TRUE,
                            doy.break = NULL, b2e.con = FALSE, t.grid.num = 50,
                            clipRange = c(0, 1800), clipMethod = "nnr"){
  if(intermediate.save){
    if(!dir.exists(intermediate.dir)){
      cat(paste0("Folder", intermediate.dir, "is created to save intermediate results."))
      dir.create(intermediate.dir, recursive = TRUE)
    }
  }
  
  imat = mat
  ###################################
  #### 1. Overall mean estimaton ####
  ###################################
  if(use.intermediate.result & file.exists(paste0(intermediate.dir, "meanest.rds"))){
    meanest = readRDS(paste0(intermediate.dir, "meanest.rds"))
  } else {
    meanest = meanEst(doy, mat, doyeval = doyeval, outlier.tol = outlier.tol, clipRange = clipRange,
                      clipMethod = clipMethod, img.nrow = img.nrow, img.ncol = img.ncol)
    if(intermediate.save)
      saveRDS(meanest, paste0(intermediate.dir, "meanest.rds"))
  }
  
  ###################################
  #### 2. Time effect estimation ####
  ###################################
  ## remove outlier pixels
  for(i in 1:length(meanest$outlier$outidx)){
    mat[meanest$outlier$outidx[i], meanest$outlier$outlst[[i]]] = NA
  }
  ## remove outlier images
  outlier.img.idx = meanest$idx$idx.outlier
  for(i in outlier.img.idx){
    mat[outlier.img.idx,] = NA
  }
  ## calculate the residuals
  rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
  
  ## estimate the temporal effect using residuals
  ## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
  if(teff){
    if(use.intermediate.result & file.exists(paste0(intermediate.dir, "teffarray.rds"))){
      teffarray = readRDS(paste0(intermediate.dir, "teffarray.rds"))
    } else {
      teffarray = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = h.tcov, h.sigma2 = h.tsigma2, t.grid.num = t.grid.num)
      if(intermediate.save)
        saveRDS(teffarray, paste0(intermediate.dir, "teffarray.rds"))
    }
    
    
    ######################################
    #### 3. Spatial effect estimation ####
    ######################################
    ## calculate residuals after removing temporal effect
    yearidx = unlist(lapply(year, function(x, y)
      which(y == x), y = as.numeric(dimnames(teffarray)[[1]])))
    doyidx = unlist(lapply(doy, function(x, y)
      which(y == x), y = as.numeric(dimnames(teffarray)[[2]])))
    for (i in 1:nrow(rmat)) {
      rmat[i, ] = rmat[i, ] - teffarray[yearidx[i], doyidx[i], ]
    }
  }
  
  ## estimate the spatial effect using residuals
  ## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
  if(seff){
    if(use.intermediate.result & file.exists(paste0(intermediate.dir, "seffmat.rds"))){
      seffmat = readRDS(paste0(intermediate.dir, "seffmat.rds"))
    } else {
      if(is.null(doy.break)){
        seffmat = seffEst(rmat, img.nrow, img.ncol, nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)$seffmat
      } else {
        seffmat = matrix(0, nrow(rmat), ncol(rmat))
        if(b2e.con){
          tmpidx = (doy <= doy.break[1]) | (doy > doy.break[length(doy.break)])
          seffmat[tmpidx,] = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                                                  nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)$seffmat
          for(i in 2:length(doy.break)){
            tmpidx = (doy <= doy.break[i]) & (doy > doy.break[i-1])
            seffmat[tmpidx,] = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                                       nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)$seffmat
          }
        } else {
          seffmat[doy <= doy.break[1],] = seffEst(rmat[doy <= doy.break[1],], img.nrow, img.ncol, 
                                                  nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)$seffmat
          for(i in 2:length(doy.break)){
            tmpidx = (doy <= doy.break[i]) & (doy > doy.break[i-1])
            seffmat[tmpidx,] = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                                                    nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)$seffmat
          }
          if(doy.break[i] < max(doy)){
            seffmat[doy > doy.break[i],] = seffEst(rmat[doy > doy.break[i],], img.nrow, img.ncol, 
                                                    nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)$seffmat
          }
        }
      }
      if(intermediate.save)
        saveRDS(seffmat, paste0(intermediate.dir, "seffmat.rds"))
    }
  }
  
  #######################
  #### 4. Gapfilling ####
  #######################
  ## partially missing images: mean + time effect + spatial effect 
  ## all missing images: mean + time effect
  ## first calculate the theoretically imputed mat.
  mat_imputed = meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
  if(teff){
    for(i in 1:nrow(mat_imputed)){
      mat_imputed[i,] = mat_imputed[i,] + teffarray[yearidx[i], doyidx[i],]
    }
  }
  if(seff){
    mat_imputed = mat_imputed + seffmat
  }
  
  #### final imputation
  outlier.action = match.arg(outlier.action)
  if(outlier.action == "keep")
    ## if keep originally observed values
    imat[is.na(imat)] = mat_imputed[is.na(imat)] else
      if (outlier.action == "remove") {
        ## if remove outliers
        imat = mat
        imat[is.na(imat)] = mat_imputed[is.na(imat)]
      }
  return(list(imat = imat, idx = meanest$idx, outlier = meanest$outlier))
}


