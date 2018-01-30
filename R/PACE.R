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
    sigma2 = max(0.5, sigma2)
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
