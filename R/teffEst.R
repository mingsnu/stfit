#' Temporal effect estimation
#'
#' @param ids ids for 'group', for data with repeated measurement over years, year is ids; 
#' for pixels belong to certain clusters, cluster is ids.
#' @param doy vector for day of year
#' @param rmat residual matrix with rows corresponding to `doy` and columns corresponding to pixel index
#' @param h.cov bandwidth for temporal covariance estimation; ignored if `weight.cov` is supplied
#' @param h.sigma2 bandwidth for sigma2 estimation
#' @param weight.cov weight vector for temporal covariance estimation
#' @param weight.sigma2 weight vector for sigma2 estimation; not used for now...
#' @param pve percentage of variance explained; used for number of eigen values selection
#' @param doyeval a vector on which to estimate
#'
#' @return
#' @export
#'
#' @examples
teffEst <- function(ids, doy, rmat,
                    doyeval = seq(min(doy), max(doy)), h.cov = 100, h.sigma2 = 300,
                     weight.cov = NULL, weight.sigma2 = NULL,
                    pve = 0.99){
  cat("Estimating the temporal effect...\n")
  if(is.null(weight.cov))
    weight.cov = weightVector(h.cov)
  if(is.null(weight.sigma2))
    weight.sigma2 = weightVector(h.sigma2)
  yeareval = sort(unique(ids))
  acomb <- function(...) abind::abind(..., along=3)
  teffarray = foreach(i = 1:ncol(rmat), .combine = "acomb", .multicombine=TRUE) %dopar% {
    resid = rmat[,i]
    nnaidx = !is.na(resid)
    if(sum(nnaidx) == 0)
      return(matrix(NA, length(yeareval), length(doyeval)))
    R0.hat = lc_cov_1d_est(ids[nnaidx], doy[nnaidx], resid[nnaidx], weight.cov, doyeval)
    # persp(tt,tt, R0.hat, theta=30, phi=30, expand=0.5, col='lightblue',
    #       xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
    sigma2 = llreg(doy, resid^2, x.eval = doyeval, h = h.sigma2)
    
    ## eigen decomposition
    ev = eigen(R0.hat)
    ev$values[ev$values < 0] = 0
    ev.idx = max(which.min(cumsum(ev$values)/sum(ev$values) < pve), 2) ## select at least 2 eigen components
    ## cat("The first ", ev.idx, " eigen values are used...\n")
    ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
    ev.val = ev$values[1:ev.idx]
    nugg = mean(sigma2 - diag(R0.hat))
    res = PACE1d(ids[nnaidx], doy[nnaidx], resid[nnaidx], ev.vec, nugg, ev.val, doyeval, yeareval)
    res$mat
  }
  dimnames(teffarray) = list(yeareval, doyeval, 1:ncol(rmat))
  return(teffarray)
}