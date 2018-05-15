#' The following function deals with the imputation using sparse pca. 
#' 
#' @param mat data matrix. Each row is a row stacked image
#' @param ev.vec eigen vectors
#' @param sigma2 measurement error
#' @param ev.val eigen values
#' 
#' @return an imputed mat matrix
PACE = function(mat, ev.vec, sigma2, ev.val){
  xi.mat = foreach (i = 1:nrow(mat), .combine = cbind) %dopar%{
    nonna.idx = !is.na(mat[i,])
    non.na.num = sum(nonna.idx)
    ## could be a problem when sigma2 is small
    ## sigma2 = max(0.5, sigma2)
    sigma2 = max(0.0001, sigma2)
    diag.mat = diag(sigma2, non.na.num, non.na.num)
    if(length(ev.val) == 1){
      xi.result = ev.val * t(ev.vec[nonna.idx, ]) %*%
        solve(ev.val* ev.vec[nonna.idx, ] %*% t(ev.vec[nonna.idx, ]) +  diag.mat) %*%
        mat[i,nonna.idx]
    } else{
      phi = ev.vec[nonna.idx,, drop=FALSE]
      tmp = t(phi) * ev.val
      xi.result =tmp %*% solve(phi %*% tmp + diag.mat, mat[i,nonna.idx])
    }
    xi.result
  }
  t(ev.vec %*% xi.mat)
}

PACE1d = function(ids, doy, resid, ev.vec, nugg, ev.val, doyeval, idseval){
  if(missing(idseval))
    idseval = sort(unique(ids))
  if(!all(idseval %in% ids))
    stop("idseval is not in ids.")
  xi.mat = foreach(i = 1:length(idseval), .combine = cbind) %dopar%{
    id.idx = which(ids == idseval[i])
    doy.idx = which(doyeval %in% doy[id.idx])
    if(length(ev.val) == 1){
      xi.result = ev.val * t(ev.vec[doy.idx, ]) %*%
        solve(ev.val* ev.vec[doy.idx, ] %*% t(ev.vec[doy.idx, ]) +  diag(nugg, length(id.idx))) %*% resid[id.idx]
    } else{
      phi = ev.vec[doy.idx,, drop=FALSE]
      tmp = t(phi) * ev.val
      xi.result =tmp %*% solve(phi %*% tmp + diag(nugg, length(id.idx)), resid[id.idx])
    }
    xi.result
  }
  return(list(idseval = idseval,
              mat = t(ev.vec %*% xi.mat)))
}

