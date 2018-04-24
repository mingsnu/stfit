#' Title
#'
#' @param doy 
#' @param mat 
#' @param img.nrow 
#' @param img.ncol 
#' @param doyeval 
#' @param h.tcov 
#' @param h.tsigma2 
#' @param h.scov 
#' @param h.ssigma2 
#' @param nnr 
#' @param msk vector of mask
#' @param ncluster number of clusters. If 0 no cluster analysis is performed and the mean is estimated based on
#' pixel level.
#' @param max.per.cluster maximum number of pixels to select for covariance estimation.
#' @param t.grid.num number of grid on doyeval, on which to estimate the covariance matrix.
#' @param breaks if NULL, estimate spatial effect on the whole image; if a list with 'block.nrow', 'block.ncol',
#' 'img.nrow' and 'img.ncol', the image is broken down into small blocks and do estimation seperately.
#' @param outlier.action 
#' @param intermediate.save TRUE or FALSE. Whether to save intermediate results.
#' @param intermediate.dir Path to save intermediate results.
#' @param teff TRUE or FALSE. Whether to estimate temporal effect.
#' @param seff TRUE or FALSE. Whether to estimate spatial effect.
#'
#' @return
#' @export
#'
#' @examples
gapfill_modis <- function(doy, mat, img.nrow, img.ncol, doyeval = 1:365, h.tcov = 100, 
                                  h.tsigma2 = 300, h.scov = 2, h.ssigma2 = 2, nnr = 10, 
                                  msk = NULL, ncluster = 50, max.per.cluster = 30, t.grid.num = 50,
                                  breaks = list(block.nrow = 10, block.ncol = 10, img.nrow = 30, img.ncol = 30),
                                  outlier.action = c("keep", "remove"),
                                  intermediate.save = TRUE, intermediate.dir = "./output/", 
                                  teff = FALSE, seff = TRUE){
  
  if(intermediate.save){
    if(!dir.exists(intermediate.dir)){
      cat("Folder 'output' is created to save intermediate results.\n")
      dir.create(intermediate.dir)
    }
  }
  imat = mat
  if(is.null(msk)){
    msk = getMask(mat)
  }
  ###################################
  #### 1. Overall mean estimaton ####
  ###################################
  meanest = meanEst(doy, mat, doyeval = 1:365, msk = msk)
  if(intermediate.save)
    saveRDS(meanest, paste0(intermediate.dir, "meanest.rds"))
  if(ncluster > 0){
    #####################################################
    #### 2. Cluster analysis based on inital results ####
    #####################################################
    cmat = t(meanest$meanmat[seq(10, 350, by = 10),!msk])
    cres = kmeans(cmat, centers = ncluster, nstart = 10, iter.max = 50)
    cluster = rep(0, ncol(mat))
    cluster[!msk] = cres$cluster
    if(intermediate.save)
      saveRDS(cluster, paste0(intermediate.dir, "cluster.rds"))
    ################################################
    #### 3. Overall mean estimaton with cluster ####
    ################################################
    ## NEVER USE smooth_spline when using clusters
    meanest_cl = meanEst(doy, mat, doyeval = doyeval, cluster = cluster, msk = msk)
    if(intermediate.save)
      saveRDS(meanest_cl, paste0(intermediate.dir, "meanest_cl.rds"))
  } else
    meanest_cl = meanest

  #########################################
  #### 2. 'temporal' effect estimation ####
  #########################################
  ## remove outlier images
  outlier.img.idx = meanest_cl$idx$idx.outlier

  for(i in outlier.img.idx){
    mat[outlier.img.idx,] = NA
  }
  ## remove outlier pixels
  for(i in 1:length(meanest_cl$outlier$outidx)){
    mat[meanest_cl$outlier$outidx[i], meanest_cl$outlier$outlst[[i]]] = NA
  }
  
  ## calculate the residuals
  rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
  
  if(teff){
    ## estimate the 'temporal' effect using residuals
    cefflist = ceffEst(doy, rmat, cluster,
                        doyeval = DOYEVAL, h.cov = h.tcov, h.sigma2 = h.tsigma2,
                       max.per.cluster = max.per.cluster, t.grid.num = t.grid.num)
    if (intermediate.save)
      saveRDS(cefflist, paste0(intermediate.dir, "cefflist.rds"))
    
    ceffmat = matrix(0, nrow(rmat), ncol(rmat))
    for (i in 1:length(cefflist)) {
      ceffmat[, cefflist[[i]]$idseval] = t(cefflist[[i]]$mat)
    }
    rmat = rmat - ceffmat
    if(intermediate.save)
      saveRDS(ceffmat, paste0(intermediate.dir, "ceffmat.rds"))
  }
  ######################################
  #### 3. Spatial effect estimation ####
  ######################################
  ## estimate the spatial effect using residuals
  if(seff){
    if(is.null(breaks)){
      seffmat = seffEst(rmat, img.nrow, img.ncol, nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2, msk = msk)$sseffmat
    } else{
      nblocks = breaks$block.nrow*breaks$block.ncol
      res3.list = foreach(n=1:nblocks) %dopar% {
        ii = floor((n-1)/breaks$block.ncol) + 1
        jj = (n-1) %% breaks$block.ncol + 1
        ## block index
        bIdx = c(t(outer(seq((ii-1)*breaks$img.nrow+1, ii*breaks$img.nrow), 
                         seq((jj-1)*breaks$img.ncol+1, jj*breaks$img.ncol),
                         FUN = function(ridx, cidx){
                           (ridx-1) * img.ncol + cidx
                         })))
        seffEst(rmat[, bIdx], breaks$img.nrow, breaks$img.ncol, nnr = nnr, 
                h.cov = h.scov, h.sigma2 = h.ssigma2, msk = msk[bIdx])
      }
      
      seffmat = matrix(NA, nrow(rmat), ncol(rmat))
      for(n in 1:nblocks){
        ii = floor((n-1)/breaks$block.ncol) + 1
        jj = (n-1) %% breaks$block.ncol + 1
        ## block index
        bIdx = c(t(outer(seq((ii-1)*breaks$img.nrow+1, ii*breaks$img.nrow), 
                         seq((jj-1)*breaks$img.ncol+1, jj*breaks$img.ncol),
                         FUN = function(ridx, cidx){
                           (ridx-1) * img.ncol + cidx
                         })))
        seffmat[,bIdx] = res3.list[[n]]$seffmat
      }
    }
    if(intermediate.save)
      saveRDS(seffmat, paste0(intermediate.dir, "seffmat.rds"))
  }

  #######################
  #### 4. Gapfilling ####
  #######################
  ## partially missing images: mean + time effect + spatial effect 
  ## all missing images: mean + time effect
  ## first calculate the theoretically imputed mat.
  mat_imputed = meanest$meanmat
  if(teff)
    mat_imputed = mat_imputed + ceffmat
  if(seff)
    mat_imputed = mat_imputed + seffmat
  
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
  return(imat)
}
