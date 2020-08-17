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
# stfit  146.3683  170.3213  182.5146  172.0806  181.5236  242.279     5