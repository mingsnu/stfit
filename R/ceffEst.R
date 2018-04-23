#' Temporal effect estimation within cluster
#'
#' @param doy vector for day of year
#' @param rmat residual matrix with rows corresponding to `doy` and columns corresponding to pixel index
#' @param h.cov bandwidth for temporal covariance estimation; ignored if `weight.cov` is supplied
#' @param h.sigma2 bandwidth for sigma2 estimation
#' @param weight.cov weight vector for temporal covariance estimation
#' @param weight.sigma2 weight vector for sigma2 estimation; not used for now...
#' @param pve percentage of variance explained; used for number of eigen values selection
#' @param doyeval a vector on which to estimate
#' @param cluster non-negative integer numbers representing the cluster for each column of rmat belongs to.
#' 0 represents mask.
#' @param max.per.cluster maximum number of pixels used tfor covariance estimation
#' @param t.grid vector on which to estimate covariance matrix
#' @param t.grid.num number of grid points on the interval of doyeval to use for covariance estimation. not used if t.grid is given.
#'
#' @return
#' @export
#'
#' @examples
ceffEst <- function(doy, rmat, cluster,
                    doyeval = seq(min(doy), max(doy)), h.cov = 100, h.sigma2 = 300,
                    weight.cov = NULL, weight.sigma2 = NULL, max.per.cluster = 30,
                    pve = 0.99, t.grid = NULL, t.grid.num = 50){
  if(length(doy) != nrow(rmat))
    stop("length of doy does not match number of row of rmat.")
  if(is.null(weight.cov))
    weight.cov = weightVector(h.cov)
  if(is.null(weight.sigma2))
    weight.sigma2 = weightVector(h.sigma2)
  cl.uni = sort(unique(cluster[cluster > 0]))
  ## t.grid: grid on which to calculate the covariance matrix
  if(is.null(t.grid)){
    if(t.grid.num < length(doyeval))
      t.grid = doyeval[unique(round(seq(1, length(doyeval), length.out=t.grid.num)))] else
        stop("t.grid.num is bigger than length of doyeval.")
  }
    
  cefflist = foreach(i = 1:length(cl.uni)) %dopar% {
    idx0 = which(cluster == cl.uni[i])
    #### sample some pixels for the covariance estimation;
    if(length(idx0) > max.per.cluster)
      idx = sample(idx0, max.per.cluster)
    resid = c(rmat[, idx])
    nnaidx = !is.na(resid)
    if(sum(nnaidx) == 0)
      return(matrix(NA, length(cl.uni), length(doyeval)))
    ids = rep(1:length(idx), each = length(doy))
    time = rep(doy, length(idx))
    
    R0.hat = lc_cov_1d_est(ids[nnaidx], time[nnaidx], resid[nnaidx], weight.cov, t.grid)
    # persp(t.grid, t.grid, R0.hat, theta=30, phi=30, expand=0.5, col='lightblue',
    #       xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
    sigma2 = llreg(time[nnaidx], resid[nnaidx]^2, x.eval = t.grid, h = h.sigma2)
    
    #### eigen decomposition
    ev = eigen(R0.hat)
    ev$values[ev$values < 0] = 0
    ev.idx = max(which.min(cumsum(ev$values)/sum(ev$values) < pve), 2) ## select at least 2 eigen components
    ## cat("The first ", ev.idx, " eigen values are used...\n")
    ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
    ev.val = ev$values[1:ev.idx]
    ## interpolate eigen function on dense grids 'doyeval'
    ev.vec = phi.interp(doyeval, ev.vec, t.grid)
    # aa = ev.vec %*% diag(ev.val)%*% t(ev.vec)
    # persp(doyeval, doyeval, aa, theta=30, phi=30, expand=0.5, col='lightblue',
    #         xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
    ## estimate the nugget effect on 1/4 ~ 3/4 time interval
    idx1 = round(length(t.grid)/4):round(length(t.grid)/4*3)
    nugg = max(mean((sigma2[idx1]-diag(R0.hat)[idx1])), 0)
    resid = c(rmat[, idx0])
    nnaidx = !is.na(resid)
    ids = rep(idx0, each = length(doy))
    time = rep(doy, length(idx0))
    res = PACE1d(ids[nnaidx], time[nnaidx], resid[nnaidx], ev.vec, nugg, ev.val, doyeval, idx0)
    res
  }
  return(cefflist)
}
