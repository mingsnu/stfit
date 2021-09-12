library(xtable)
printbold <-
  function(x, which = NULL, each = c("column", "row"), max = TRUE,
           NA.string = "", type = c("latex", "html"),
           sanitize.text.function = force,
           sanitize.rownames.function = NULL,
           sanitize.colnames.function = NULL, ...)
  {
    stopifnot(inherits(x, "xtable"))
    each <- match.arg(each)
    type <- match.arg(type)
    digits <- rep(digits(x), length = ncol(x)+1)
    if (!is.null(which)) {
      stopifnot(nrow(which) == nrow(x))
      stopifnot(ncol(which) == ncol(x))
      boldmatrix <- which
    } else {
      boldmatrix <- matrix(FALSE, ncol = ncol(x), nrow = nrow(x))
      ## round values before calculating max/min to avoid trivial diffs
      for (i in 1:ncol(x)) {
        if (!is.numeric(x[,i])) next
        x[,i] <- round(x[,i], digits = digits[i+1])
      }
      if (each == "column") {
        max <- rep(max, length = ncol(x))
        for (i in 1:ncol(x)) {
          xi <- x[,i]
          if (!is.numeric(xi)) next
          if (is.na(max[i])) next
          imax <- max(xi, na.rm = TRUE)
          if (!max[i])
            imax <- min(xi, na.rm = TRUE)
          boldmatrix[xi == imax, i] <- TRUE
        }
      } else if (each == "row") {
        max <- rep(max, length = nrow(x))
        for (i in 1:nrow(x)) {
          xi <- x[i,]
          ok <- sapply(xi, is.numeric)
          if (!any(ok)) next
          if (is.na(max[i])) next
          imax <- max(unlist(xi[ok]), na.rm = TRUE)
          if (!max[i])
            imax <- min(unlist(xi[ok]), na.rm = TRUE)
          whichmax <- sapply(xi, identical, imax)
          boldmatrix[i, whichmax] <- TRUE
        }
      }
    }
    ## need to convert to character
    ## only support per-column formats, not cell formats
    display <- rep(display(x), length = ncol(x)+1)
    for (i in 1:ncol(x)) {
      if (!is.numeric(x[,i])) next
      ina <- is.na(x[,i])
      x[,i] <- formatC(x[,i], digits = digits[i+1],
                       format = display[i+1])
      x[ina, i] <- NA.string
      display(x)[i+1] <- "s"
      ## embolden
      yes <- boldmatrix[,i]
      if (type == "latex") {
        x[yes,i] <- paste("\\textbf{", x[yes,i], "}", sep = "")
      } else {
        x[yes,i] <- paste("<strong>", x[yes,i], "</strong>", sep = "")
      }
    }
    print(x, ..., type = type, NA.string = NA.string,
          sanitize.text.function = sanitize.text.function,
          sanitize.rownames.function = sanitize.rownames.function,
          sanitize.colnames.function = sanitize.colnames.function)
  }


####################################################
############ Load all experiment results ###########
####################################################
stfit_rmse_A = do.call("cbind", 
                     lapply(paste0("./landsat_simulation_study/stfit_siteA/", 
                                   c("stfit_spring", "stfit_summer", 
                                     "stfit_fall", "stfit_winter"), "/res.rds"),
                            function(f) {
                              matrix(unlist(lapply(readRDS(f), function(x)x[1])), 15)
                            }))
gapfill_rmse_A = do.call("cbind", 
                       lapply(paste0("./landsat_simulation_study/gapfill_siteA/", 
                                     c("gapfill_spring", "gapfill_summer", 
                                       "gapfill_fall", "gapfill_winter"), "/res.rds"),
                          function(f) {
                            matrix(unlist(lapply(readRDS(f), function(x)x[1])), 15)
                          }))
krig_detrended_rmse_A = readRDS("./landsat_simulation_study/kriging_detrended_siteA/output/res.rds")[[1]]
krig_not_detrended_rmse_A = readRDS("./landsat_simulation_study/kriging_not_detrended_siteA/output/res.rds")[[1]]

stfit_rmse_B = do.call("cbind", 
                       lapply(paste0("./landsat_simulation_study/stfit_siteB/", 
                                     c("stfit_spring", "stfit_summer", 
                                       "stfit_fall", "stfit_winter"), "/res.rds"),
                              function(f) {
                                matrix(unlist(lapply(readRDS(f), function(x)x[1])), 15)
                              }))
gapfill_rmse_B = do.call("cbind", 
                         lapply(paste0("./landsat_simulation_study/gapfill_siteB/", 
                                       c("gapfill_spring", "gapfill_summer", 
                                         "gapfill_fall", "gapfill_winter"), "/res.rds"),
                                function(f) {
                                  matrix(unlist(lapply(readRDS(f), function(x)x[1])), 15)
                                }))
krig_detrended_rmse_B = readRDS("./landsat_simulation_study/kriging_detrended_siteB/output/res.rds")[[1]]
krig_not_detrended_rmse_B = readRDS("./landsat_simulation_study/kriging_not_detrended_siteB/output/res.rds")[[1]]

########################################################
###### Table S3: Three methods' RMSEs for spring #######
########################################################
RMSEmat1 = stfit_rmse_A[,1:5]
RMSEmat2 = gapfill_rmse_A[,1:5]
RMSEmat3 = krig_not_detrended_rmse_A[,1:5]
RMSEmat4 = krig_detrended_rmse_A[,1:5]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S3.txt")

########################################################
###### Table S4: Three methods' RMSEs for summer #######
########################################################
RMSEmat1 = stfit_rmse_A[,6:10]
RMSEmat2 = gapfill_rmse_A[,6:10]
RMSEmat3 = krig_not_detrended_rmse_A[,6:10]
RMSEmat4 = krig_detrended_rmse_A[,6:10]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S4.txt")

######################################################
###### Table S5: Three methods' RMSEs for fall #######
######################################################
RMSEmat1 = stfit_rmse_A[,11:15]
RMSEmat2 = gapfill_rmse_A[,11:15]
RMSEmat3 = krig_not_detrended_rmse_A[,11:15]
RMSEmat4 = krig_detrended_rmse_A[,11:15]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S5.txt")

########################################################
###### Table S6: Three methods' RMSEs for winter #######
########################################################
RMSEmat1 = stfit_rmse_A[,16:19]
RMSEmat2 = gapfill_rmse_A[,16:19]
RMSEmat3 = krig_not_detrended_rmse_A[,16:19]
RMSEmat4 = krig_detrended_rmse_A[,16:19]
RMSEmat = matrix(0, 15, 16)
RMSEmat[, seq(1, 16, 4)] = RMSEmat1
RMSEmat[, seq(2, 16, 4)] = RMSEmat2
RMSEmat[, seq(3, 16, 4)] = RMSEmat3
RMSEmat[, seq(4, 16, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 16)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 16, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 16, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 16, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 16, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S6.txt")

########################################################
###### Table S7: Three methods' RMSEs for spring #######
########################################################
RMSEmat1 = stfit_rmse_B[,1:5]
RMSEmat2 = gapfill_rmse_B[,1:5]
RMSEmat3 = krig_not_detrended_rmse_B[,1:5]
RMSEmat4 = krig_detrended_rmse_B[,1:5]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S7.txt")

########################################################
###### Table S8: Three methods' RMSEs for summer #######
########################################################
RMSEmat1 = stfit_rmse_B[,6:10]
RMSEmat2 = gapfill_rmse_B[,6:10]
RMSEmat3 = krig_not_detrended_rmse_B[,6:10]
RMSEmat4 = krig_detrended_rmse_B[,6:10]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S8.txt")

######################################################
###### Table S9: Three methods' RMSEs for fall #######
######################################################
RMSEmat1 = stfit_rmse_B[,11:15]
RMSEmat2 = gapfill_rmse_B[,11:15]
RMSEmat3 = krig_not_detrended_rmse_B[,11:15]
RMSEmat4 = krig_detrended_rmse_B[,11:15]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S9.txt")


########################################################
###### Table S10: Three methods' RMSEs for winter #######
########################################################
RMSEmat1 = stfit_rmse_B[,16:20]
RMSEmat2 = gapfill_rmse_B[,16:20]
RMSEmat3 = krig_not_detrended_rmse_B[,16:20]
RMSEmat4 = krig_detrended_rmse_B[,16:20]
RMSEmat = matrix(0, 15, 20)
RMSEmat[, seq(1, 20, 4)] = RMSEmat1
RMSEmat[, seq(2, 20, 4)] = RMSEmat2
RMSEmat[, seq(3, 20, 4)] = RMSEmat3
RMSEmat[, seq(4, 20, 4)] = RMSEmat4
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = rep(c('stfit', 'gapfill', 'kriging', 'kriging (detrended)'),5)
RMSEmat
## for latex output
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 20)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 20, 4)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3) & (RMSEmat1 < RMSEmat4)
boldmat[, seq(2, 20, 4)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3) & (RMSEmat2 < RMSEmat4)
boldmat[, seq(3, 20, 4)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2) & (RMSEmat3 < RMSEmat4)
boldmat[, seq(4, 20, 4)] = (RMSEmat4 < RMSEmat1) & (RMSEmat4 < RMSEmat2) & (RMSEmat4 < RMSEmat3)
printbold(RMSExtable, which = boldmat, NA.string = "-",
          include.rownames=TRUE, file = "output/tab_S10.txt")

###################################################################
###### Table 1: Average RMSE by missing pattern and seasons #######
###################################################################
dfB = landsat106 %>% filter(year >= 2000)
matB = as.matrix(dfB[,-c(1:2)])
pidx = c(68, 209, 352, 605, 624, 74, 156, 263, 273, 499, 184, 369, 508, 517, 565)
pmat = matB[pidx,]
pmat[!is.na(pmat)] = 1
n.na = apply(pmat, 1, function(x) sum(is.na(x)))
n.na /961 * 100

## site A
stfit_rmse_int_A = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    if(j < 4){
      tmpvec = c(stfit_rmse_A[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
      weight = rep(n.na[seq((i-1)*5+1, i*5)], 5)
      nonna.idx = !is.na(tmpvec)
      stfit_rmse_int_A[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
    } else {
      tmpvec = c(stfit_rmse_A[seq((i-1)*5+1, i*5), 16:19])
      weight = rep(n.na[seq((i-1)*5+1, i*5)], 4)
      nonna.idx = !is.na(tmpvec)
      stfit_rmse_int_A[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
    }
  }
}
rownames(stfit_rmse_int_A) = c("< 30%", "30%~70%", "70%~99%")
colnames(stfit_rmse_int_A) = c("Spring", "Summer", "Fall", "Winter")

## site B
stfit_rmse_int_B = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
      tmpvec = c(stfit_rmse_B[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
      weight = rep(n.na[seq((i-1)*5+1, i*5)], 5)
      nonna.idx = !is.na(tmpvec)
      stfit_rmse_int_B[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
  }
}
rownames(stfit_rmse_int_B) = c("< 30%", "30%~70%", "70%~99%")
colnames(stfit_rmse_int_B) = c("Spring", "Summer", "Fall", "Winter")

cbind(stfit_rmse_int_A, stfit_rmse_int_B)
print(xtable(cbind(stfit_rmse_int_A, stfit_rmse_int_B)),
      include.rownames=FALSE, file = "output/tab_1.txt")

#####################################################################
###### Table 3 (a): Average RE by missing pattern and seasons #######
#####################################################################
RE_21 = gapfill_rmse_A/stfit_rmse_A
RE_31 = krig_not_detrended_rmse_A/stfit_rmse_A
RE_41 = krig_detrended_rmse_A/stfit_rmse_A
#### RE matrix for gapfill
RE2 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    if(j < 4){
      RE2[i,j] = mean(RE_21[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm = TRUE)
    } else {
      RE2[i,j] = mean(RE_21[seq((i-1)*5+1, i*5), 16:19], na.rm = TRUE)
    }
  }
}

#### RE matrix for kriging not detrended
RE3 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    if(j < 4){
      RE3[i,j] = mean(RE_31[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm = TRUE)
    } else {
      RE3[i,j] = mean(RE_31[seq((i-1)*5+1, i*5), 16:19], na.rm = TRUE)
    }
  }
}

#### RE matrix for kriging detrended
RE4 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    if(j < 4){
      RE4[i,j] = mean(RE_41[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm = TRUE)
    } else {
      RE4[i,j] = mean(RE_41[seq((i-1)*5+1, i*5), 16:19], na.rm = TRUE)
    }
  }
}

RE = matrix(NA, 3, 12)
RE[, seq(1,12, by=3)] = RE2
RE[, seq(2,12, by=3)] = RE3
RE[, seq(3,12, by=3)] = RE4
rownames(RE) = c("< 30%", "30%~70%", "70%~99%")
colnames(RE) = rep(c("Spring", "Summer", "Fall", "Winter"), each=3)
RE
print(xtable(RE), include.rownames=FALSE, file = "output/tab_3a.txt")

#####################################################################
###### Table 3 (b): Average RE by missing pattern and seasons #######
#####################################################################
RE_21 = gapfill_rmse_B/stfit_rmse_B
RE_31 = krig_not_detrended_rmse_B/stfit_rmse_B
RE_41 = krig_detrended_rmse_B/stfit_rmse_B
#### RE matrix for gapfill
RE2 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
      RE2[i,j] = mean(RE_21[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm = TRUE)
  }
}
#### RE matrix for kriging not detrended
RE3 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
      RE3[i,j] = mean(RE_31[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm = TRUE)
  }
}
#### RE matrix for kriging detrended
RE4 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    RE4[i,j] = mean(RE_41[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm = TRUE)
  }
}

RE = matrix(NA, 3, 12)
RE[, seq(1,12, by=3)] = RE2
RE[, seq(2,12, by=3)] = RE3
RE[, seq(3,12, by=3)] = RE4
rownames(RE) = c("< 30%", "30%~70%", "70%~99%")
colnames(RE) = rep(c("Spring", "Summer", "Fall", "Winter"), each=3)
RE
print(xtable(RE), include.rownames=FALSE, file = "output/tab_3b.txt")

###############################################################
###### Table 2: Contribution of temporal/spatial effect #######
###############################################################
rmse_mean = matrix(unlist(
  lapply(readRDS("./landsat_simulation_study/effects_siteA/output/mean_res.rds"), function(x)x[1])), 15)
rmse_teff = matrix(unlist(
  lapply(readRDS("./landsat_simulation_study/effects_siteA/output/teff_res.rds"), function(x)x[1])), 15)
rmse_seff = matrix(unlist(
  lapply(readRDS("./landsat_simulation_study/effects_siteA/output/seff_res.rds"), function(x)x[1])), 15)
rmse_steff = do.call("cbind", 
                     lapply(paste0("./landsat_simulation_study/stfit_siteA/", 
                                   c("stfit_spring", "stfit_summer", 
                                     "stfit_fall", "stfit_winter"), "/res.rds"),
                            function(f) {
                              matrix(unlist(lapply(readRDS(f), function(x)x[1])), 15)
                            }))

rrmse_t = rmse_teff/rmse_mean
rrmse_s = rmse_seff/rmse_mean
rrmse_st = rmse_steff/rmse_mean
rrmse_A = c(mean(rrmse_t), mean(rrmse_s), mean(rrmse_st))
rrmse_A

rmse_mean = matrix(unlist(
  lapply(readRDS("./landsat_simulation_study/effects_siteB/output/mean_res.rds"), function(x)x[1])), 15)
rmse_teff = matrix(unlist(
  lapply(readRDS("./landsat_simulation_study/effects_siteB/output/teff_res.rds"), function(x)x[1])), 15)
rmse_seff = matrix(unlist(
  lapply(readRDS("./landsat_simulation_study/effects_siteB/output/seff_res.rds"), function(x)x[1])), 15)
rmse_steff = do.call("cbind", 
                     lapply(paste0("./landsat_simulation_study/stfit_siteB/", 
                                   c("stfit_spring", "stfit_summer", 
                                     "stfit_fall", "stfit_winter"), "/res.rds"),
                            function(f) {
                              matrix(unlist(lapply(readRDS(f), function(x)x[1])), 15)
                            }))
rrmse_t = rmse_teff/rmse_mean
rrmse_s = rmse_seff/rmse_mean
rrmse_st = rmse_steff/rmse_mean
rrmse_B = c(mean(rrmse_t), mean(rrmse_s), mean(rrmse_st))
rrmse_B

eff_mat = rbind(rrmse_A, rrmse_B)
rownames(eff_mat) = c("A", "B")
colnames(eff_mat) = c("T", "S", "TS")
eff_mat
print(xtable(eff_mat),include.rownames=FALSE, file = "output/tab_2.txt")
