library(microbenchmark)
library(stfit)
library(dplyr)
library(doParallel)

df = landsat2 %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])

doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))
datarray = array(NA, dim = c(31, 31, 46, 16), dimnames = list(1:31, 1:31, doybinuni, yearuni))
for(ii in 1:16){
  for(jj in 1:46){
    idx = year == yearuni[ii] & doybin == doybinuni[jj]
    if(sum(idx) == 1)
      datarray[,,jj,ii] = matrix(mat[year == yearuni[ii] & doybin == doybinuni[jj],], 31) else
        if(sum(idx) > 1)
          warning("Multiple matches.")
  }
}

.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X[x.eval,] %*% lmfit$coefficient)
}
stfit::opts_stfit$set(temporal_mean_est = customfun)
registerDoParallel(16)
microbenchmark("gapfill" = {
  res1 = gapfill::Gapfill(datarray, clipRange = c(0, 3000), dopar = TRUE)
},
"stfit" = {
  res2 <- stfit_landsat(year, doy, mat, 31, 31, nnr=30, clipRange= c(0,3000),
                        use.intermediate.result = FALSE, intermediate.save = FALSE)
},
times = 5)
# Unit: seconds
# expr       min        lq      mean    median        uq      max neval
# gapfill 1552.5008 1580.4793 1769.5675 1633.1508 1645.4840 2436.223     5
# stfit   115.267   115.4007  117.8956  117.7408  118.3663  122.7031     5