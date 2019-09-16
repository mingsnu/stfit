# customized function for estimating pooled temporal and cluster effects

cteffEst <- function(year, doy, rmat, cluster,
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
  yeareval <- sort(unique(year))
  ## t.grid: grid on which to calculate the covariance matrix
  if(is.null(t.grid)){
    if(t.grid.num < length(doyeval))
      t.grid = doyeval[unique(round(seq(1, length(doyeval), length.out=t.grid.num)))] else
        stop("t.grid.num is bigger than length of doyeval.")
  }
  
  ctefflist = foreach(i = 1:length(cl.uni)) %dopar% {
    # pixels of cluster i
    idx0 = which(cluster == cl.uni[i])
    #### sample some pixels for the covariance estimation;
    if(length(idx0) > max.per.cluster)
      idx = sample(idx0, max.per.cluster)
    # residual matrix: row is pixels, transformed to vector
    resid = c(rmat[, idx])
    # nonmissing residuals
    nnaidx = !is.na(resid)
    
    if(sum(nnaidx) == 0) {
      res <- matrix(NA, length(idx0), length(doy))
      return(list(mat = res, pixels = idx0))
    }
    ids <- rep(year, length(idx))
    #ids = rep(1:length(idx), each = length(doy))
    time = rep(doy, length(idx))
    
    # covariance function at t.grid
    R0.hat = lc_cov_1d_est(ids[nnaidx], time[nnaidx], resid[nnaidx], weight.cov, t.grid)
    # persp(t.grid, t.grid, R0.hat, theta=30, phi=30, expand=0.5, col='lightblue',
    #       xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
    sigma2 = llreg(time[nnaidx], resid[nnaidx]^2, x.eval = t.grid, h = h.sigma2)
    
    #### eigen decomposition
    cat("Doing eigen decomposition.\n")
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
    res <- matrix(0, nrow = length(idx0), ncol = length(doy))
    # get fpca result for each pixel within the current cluster
    cat("Calculating imputation results.\n")
    for (i in 1:length(idx0)) {
      resid <- rmat[, idx0[i]]
      nnaidx <- !is.na(resid)
      if(!all(yeareval %in% year[nnaidx])){
        # estimated effect with yeareval as row, and doyeval as column
        resi = PACE1d(year[nnaidx], doy[nnaidx], resid[nnaidx], ev.vec, nugg, ev.val, doyeval)
        eff.mat <- resi$mat
        yrevals <- resi$idseval
        eff.list <- lapply(1:nrow(eff.mat), function(x) {eff.mat[x, which(doyeval %in% doy[year == yrevals[x]])]})
        res[i, which(year %in% yrevals)] <- unlist(eff.list)
      } else{
        resi = PACE1d(year[nnaidx], doy[nnaidx], resid[nnaidx], ev.vec, nugg, ev.val, doyeval, yeareval)
        eff.mat <- resi$mat
        eff.list <- lapply(1:nrow(eff.mat), function(x) {eff.mat[x, which(doyeval %in% doy[year == yeareval[x]])]})
        res[i,] <- unlist(eff.list)
      }
    }
    return(list(mat = res, pixels = idx0))
  }
  return(cefflist)
}
