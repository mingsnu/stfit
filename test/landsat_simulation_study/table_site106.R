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


##########################################################################################
##########################################################################################

##################
#### including kriging
##################
RMSE_krig = readRDS("./site106_kriging/output/res.rds")[[1]]
## Spring table ===================================
res = readRDS("./site106_stfit/stfit_spring/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 15)
rownames(RMSEmat1) = paste0("P", 1:15)
colnames(RMSEmat1) = paste0("F", 1:5)
res = readRDS("./site106_gapfill/gapfill_spring/res.rds")
RMSEmat2 = matrix(unlist(lapply(res, function(x)x[1])), 15)
RMSEmat3 = RMSE_krig[,1:5]

## combine
RMSEmat = matrix(0, 15, 15)
RMSEmat[, seq(1, 15, 3)] = RMSEmat1
RMSEmat[, seq(2, 15, 3)] = RMSEmat2
RMSEmat[, seq(3, 15, 3)] = RMSEmat3
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(1:5, each=3))
RMSExtable = xtable::xtable(RMSEmat)

boldmat = matrix(TRUE, 15, 15)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 15, 3)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3)
boldmat[, seq(2, 15, 3)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3)
boldmat[, seq(3, 15, 3)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2)
printbold(RMSExtable, which = boldmat, NA.string = "-")

## Summer table ===================================
res = readRDS("./site106_stfit/stfit_summer/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 15)
rownames(RMSEmat1) = paste0("P", 1:15)
colnames(RMSEmat1) = paste0("F", 6:10)
res = readRDS("./site106_gapfill/gapfill_summer/res.rds")
RMSEmat2 = matrix(unlist(lapply(res, function(x)x[1])), 15)
RMSEmat3 = RMSE_krig[,6:10]

## combine
RMSEmat = matrix(0, 15, 15)
RMSEmat[, seq(1, 15, 3)] = RMSEmat1
RMSEmat[, seq(2, 15, 3)] = RMSEmat2
RMSEmat[, seq(3, 15, 3)] = RMSEmat3
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(6:10, each=3))
RMSExtable = xtable::xtable(RMSEmat)

boldmat = matrix(TRUE, 15, 15)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 15, 3)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3)
boldmat[, seq(2, 15, 3)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3)
boldmat[, seq(3, 15, 3)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2)
printbold(RMSExtable, which = boldmat, NA.string = "-")

## Fall table ===================================
res = readRDS("./site106_stfit/stfit_fall/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 15)
rownames(RMSEmat1) = paste0("P", 1:15)
colnames(RMSEmat1) = paste0("F", 11:15)
res = readRDS("./site106_gapfill/gapfill_fall/res.rds")
RMSEmat2 = matrix(unlist(lapply(res, function(x)x[1])), 15)
RMSEmat3 = RMSE_krig[,11:15]

## combine
RMSEmat = matrix(0, 15, 15)
RMSEmat[, seq(1, 15, 3)] = RMSEmat1
RMSEmat[, seq(2, 15, 3)] = RMSEmat2
RMSEmat[, seq(3, 15, 3)] = RMSEmat3
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(11:15, each=3))
RMSExtable = xtable::xtable(RMSEmat)

boldmat = matrix(TRUE, 15, 15)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 15, 3)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3)
boldmat[, seq(2, 15, 3)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3)
boldmat[, seq(3, 15, 3)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2)
printbold(RMSExtable, which = boldmat, NA.string = "-")

## Winter table ===================================
res = readRDS("./site106_stfit/stfit_winter/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 15)
rownames(RMSEmat1) = paste0("P", 1:15)
colnames(RMSEmat1) = paste0("F", 16:20)
res = readRDS("./site106_gapfill/gapfill_winter/res.rds")
RMSEmat2 = matrix(unlist(lapply(res, function(x)x[1])), 15)
RMSEmat3 = RMSE_krig[,16:20]

## combine
RMSEmat = matrix(0, 15, 15)
RMSEmat[, seq(1, 15, 3)] = RMSEmat1
RMSEmat[, seq(2, 15, 3)] = RMSEmat2
RMSEmat[, seq(3, 15, 3)] = RMSEmat3
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(16:20, each=3))
RMSExtable = xtable::xtable(RMSEmat)

boldmat = matrix(TRUE, 15, 15)
RMSEmat2[is.na(RMSEmat2)] = 999999
boldmat[, seq(1, 15, 3)] = (RMSEmat1 < RMSEmat2) & (RMSEmat1 < RMSEmat3)
boldmat[, seq(2, 15, 3)] = (RMSEmat2 < RMSEmat1) & (RMSEmat2 < RMSEmat3)
boldmat[, seq(3, 15, 3)] = (RMSEmat3 < RMSEmat1) & (RMSEmat3 < RMSEmat2)
printbold(RMSExtable, which = boldmat, NA.string = "-")

############################################
###### PCT stfit better than gapfill #######
############################################
res = readRDS("./site106_stfit/stfit_spring/res.rds")
RMSEmat1_spring = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_gapfill/gapfill_spring/res.rds")
RMSEmat2_spring = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_stfit/stfit_summer/res.rds")
RMSEmat1_summer = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_gapfill/gapfill_summer/res.rds")
RMSEmat2_summer = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_stfit/stfit_fall/res.rds")
RMSEmat1_fall = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_gapfill/gapfill_fall/res.rds")
RMSEmat2_fall = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_stfit/stfit_winter/res.rds")
RMSEmat1_winter = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_gapfill/gapfill_winter/res.rds")
RMSEmat2_winter = matrix(unlist(lapply(res, function(x)x[1])), 15)
RMSEmat1 = cbind(RMSEmat1_spring, RMSEmat1_summer, RMSEmat1_fall, RMSEmat1_winter)
RMSEmat2 = cbind(RMSEmat2_spring, RMSEmat2_summer, RMSEmat2_fall, RMSEmat2_winter)
RMSEmat3 = readRDS("./site106_kriging/output/res.rds")[[1]]
mat12 = RMSEmat1 < RMSEmat2
pctmat = matrix(NA, 3,4)
for(i in 1:3){
  for(j in 1:4){
    tmpvec = c(mat12[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    pctmat[i,j] = sum(tmpvec, na.rm=TRUE)/sum(!is.na(tmpvec))
  }
}
rownames(pctmat) = c("< 30%", "30%~70%", "70%~99%")
colnames(pctmat) = c("Spring", "Summer", "Fall", "Winter")
xtable::xtable(pctmat)

##########################################################
###### Average RMSE by missing pattern and seasons #######
##########################################################
df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 2000] = NA

#### partial missing image indexes with different missing percentage
pmat = readRDS("missing_pattern/output/missing_pattern.rds")
n.na = apply(pmat, 1, function(x) sum(is.na(x)))
n.na /961 * 100

#### integrated RMSE matrix for stfit
RMSEmat1.int = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    tmpvec = c(RMSEmat1[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    weight = rep(n.na[seq((i-1)*5+1, i*5)], 5)
    nonna.idx = !is.na(tmpvec)
    RMSEmat1.int[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
  }
}
RMSEmat1.int
rownames(RMSEmat1.int) = c("< 30%", "30%~70%", "70%~99%")
colnames(RMSEmat1.int) = c("Spring", "Summer", "Fall", "Winter")
xtable(RMSEmat1.int)

#### integrated RMSE matrix for gapfill
RMSEmat2.int = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    tmpvec = c(RMSEmat2[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    weight = rep(n.na[seq((i-1)*5+1, i*5)], 5)
    nonna.idx = !is.na(tmpvec)
    RMSEmat2.int[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
  }
}
RMSEmat2.int

#### integrated RMSE matrix for kriging
RMSEmat3.int = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
    tmpvec = c(RMSEmat3[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    weight = rep(n.na[seq((i-1)*5+1, i*5)], 5)
    nonna.idx = !is.na(tmpvec)
    RMSEmat3.int[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
  }
}
RMSEmat3.int

RMSEmat.int = matrix(NA, 3, 12)
RMSEmat.int[, seq(1,12, by=3)] = RMSEmat1.int
RMSEmat.int[, seq(2,12, by=3)] = RMSEmat2.int
RMSEmat.int[, seq(3,12, by=3)] = RMSEmat3.int
RMSEmat.int
rownames(RMSEmat.int) = c("< 30%", "30%~70%", "70%~99%")
colnames(RMSEmat.int) = rep(c("Spring", "Summer", "Fall", "Winter"), each=3)
RMSEmat.int
boldmat = matrix(TRUE, 3, 12)
boldmat[, seq(1,12, by=3)] = (RMSEmat1.int < RMSEmat2.int) & (RMSEmat1.int < RMSEmat3.int)
boldmat[, seq(2,12, by=3)] = (RMSEmat2.int < RMSEmat1.int) & (RMSEmat2.int < RMSEmat3.int)
boldmat[, seq(3,12, by=3)] = (RMSEmat3.int < RMSEmat1.int) & (RMSEmat3.int < RMSEmat2.int)
printbold(xtable::xtable(RMSEmat.int), which = boldmat, NA.string = "-")

##########################################################
###### Average RE by missing pattern and seasons #######
##########################################################
RE_21 = RMSEmat2/RMSEmat1
RE_31 = RMSEmat3/RMSEmat1
#### RE matrix for gapfill
RE2 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
      RE2[i,j] = mean(RE_21[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm=TRUE)
  }
}
RE2
#### RE matrix for kriging
RE3 = matrix(NA, 3, 4)
for(i in 1:3){
  for(j in 1:4){
      RE3[i,j] = mean(RE_31[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)], na.rm=TRUE)
  }
}
RE3

RE = matrix(NA, 3, 12)
RE[, seq(1,12, by=3)] = 1
RE[, seq(2,12, by=3)] = RE2
RE[, seq(3,12, by=3)] = RE3
RE
rownames(RE) = c("< 30%", "30%~70%", "70%~99%")
colnames(RE) = rep(c("Spring", "Summer", "Fall", "Winter"), each=3)
RE
boldmat = matrix(TRUE, 3, 12)
boldmat[, seq(1,12, by=3)] = (1 < RE2) & (1 < RE3)
boldmat[, seq(2,12, by=3)] = (RE2 < 1) & (RE2 < RE3)
boldmat[, seq(3,12, by=3)] = (RE3 < 1) & (RE3 < RE2)
printbold(xtable::xtable(RE), which = boldmat, NA.string = "-")

##################
##################
##### effects tables =======================
## Spring table ===================================
res = readRDS("./site106_effects/output/mean_res.rds")
RMSEmat_mean1 = matrix(unlist(lapply(res, function(x)x[1])), 15)
rownames(RMSEmat_mean1) = paste0("P", 1:15)
res = readRDS("./site106_effects/output/teff_res.rds")
RMSEmat_teff1 = matrix(unlist(lapply(res, function(x)x[1])), 15)
res = readRDS("./site106_effects/output/seff_res.rds")
RMSEmat_seff1 = matrix(unlist(lapply(res, function(x)x[1])), 15)

# 
# ## combine
# RMSEmat = matrix(0, 15, 15)
# RMSEmat[, seq(1, 15, 3)] = RMSEmat_mean1
# RMSEmat[, seq(2, 15, 3)] = RMSEmat_teff1
# RMSEmat[, seq(3, 15, 3)] = RMSEmat1
# RMSEmat
# rownames(RMSEmat) = paste0("P", 1:15)
# colnames(RMSEmat) = paste0("F", rep(1:5, each=3))
# RMSExtable = xtable::xtable(RMSEmat)
# 
# boldmat = matrix(TRUE, 15, 15)
# ## table(c(RMSEmat_mean1 < RMSEmat2))
# boldmat[, seq(1, 15, 3)] = (RMSEmat_mean1 < RMSEmat2) & (RMSEmat_mean1 < RMSEmat3)
# boldmat[, seq(2, 15, 3)] = (RMSEmat2 < RMSEmat_mean1) & (RMSEmat2 < RMSEmat3)
# boldmat[, seq(3, 15, 3)] = (RMSEmat3 < RMSEmat_mean1) & (RMSEmat3 < RMSEmat2)
# printbold(RMSExtable, which = boldmat, NA.string = "-")
RMSE_t2m_ratio = RMSEmat_teff1/RMSEmat_mean1
RMSE_s2m_ratio = RMSEmat_seff1/RMSEmat_mean1
RMSE_st2m_ratio = RMSEmat1/RMSEmat_mean1
eff_mat = cbind(apply(RMSE_t2m_ratio,1,mean), apply(RMSE_s2m_ratio,1,mean), apply(RMSE_st2m_ratio,1,mean))
eff_mat = apply(eff_mat, 2, FUN = function(x) aggregate(x,by=list(rep(1:3,each=5)),mean)$x)
eff_mat = cbind(1, eff_mat)
rownames(eff_mat) = c("< 30%", "30%~70%", "70%~99%")
colnames(eff_mat) = c("M", "MT", "MS", "MTS")
eff_mat
xtable(eff_mat)
