#### print in bold font for xtable
#### https://gist.github.com/floybix/452201
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


res = readRDS("../simulation2_nnr30/pidx0.1/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
RMSEmat2 = readRDS("./pidx0.1/RMSEmat2.rds")
res = readRDS("../simulation2_nnr30/pidx0.4_0.6/res.rds")
tmp1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
tmp2 = readRDS("./pidx0.4_0.6/RMSEmat2.rds")
RMSEmat1 = cbind(RMSEmat1, tmp1)
RMSEmat2 = cbind(RMSEmat2, tmp2)
res = readRDS("../simulation2_nnr30/pidx0.8_0.95/res.rds")
tmp1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
tmp2 = readRDS("./pidx0.8_0.95/RMSEmat2.rds")
RMSEmat1 = cbind(RMSEmat1, tmp1)
RMSEmat2 = cbind(RMSEmat2, tmp2)

rownames(RMSEmat1) = paste0("F", 1:20)
colnames(RMSEmat1) = paste0("P", 1:15)
xtable::xtable(RMSEmat1)

rownames(RMSEmat2) = paste0("F", 1:20)
colnames(RMSEmat2) = paste0("P", 1:15)
xtable::xtable(RMSEmat2)



####### combined tales
## combined table for season 1
RMSEmat = matrix(0, 15, 10)
RMSEmat[, seq(1,10, 2)] = t(RMSEmat1[1:5,])
RMSEmat[, seq(2,10, 2)] = t(RMSEmat2[1:5,])
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(1:5, each=2))
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 10)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 10, 2)] = t(RMSEmat1[1:5,] < RMSEmat2[1:5,])
boldmat[, seq(2, 10, 2)] = t(RMSEmat1[1:5,] > RMSEmat2[1:5,])
boldmat[is.na(boldmat)] = FALSE
printbold(RMSExtable, which = boldmat, NA.string = "-")

## combined table for season 2
RMSEmat = matrix(0, 15, 10)
RMSEmat[, seq(1,10, 2)] = t(RMSEmat1[6:10,])
RMSEmat[, seq(2,10, 2)] = t(RMSEmat2[6:10,])
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(6:10, each=2))
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 10)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 10, 2)] = t(RMSEmat1[6:10,] < RMSEmat2[6:10,])
boldmat[, seq(2, 10, 2)] = t(RMSEmat1[6:10,] > RMSEmat2[6:10,])
boldmat[is.na(boldmat)] = FALSE
printbold(RMSExtable, which = boldmat, NA.string = "-")

## combined table for season 3
RMSEmat = matrix(0, 15, 10)
RMSEmat[, seq(1,10, 2)] = t(RMSEmat1[11:15,])
RMSEmat[, seq(2,10, 2)] = t(RMSEmat2[11:15,])
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(11:15, each=2))
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 10)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 10, 2)] = t(RMSEmat1[11:15,] < RMSEmat2[11:15,])
boldmat[, seq(2, 10, 2)] = t(RMSEmat1[11:15,] > RMSEmat2[11:15,])
boldmat[is.na(boldmat)] = FALSE
printbold(RMSExtable, which = boldmat, NA.string = "-")

## combined table for season 4
RMSEmat = matrix(0, 15, 10)
RMSEmat[, seq(1,10, 2)] = t(RMSEmat1[16:20,])
RMSEmat[, seq(2,10, 2)] = t(RMSEmat2[16:20,])
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(16:20, each=2))
RMSExtable = xtable::xtable(RMSEmat)
boldmat = matrix(TRUE, 15, 10)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 10, 2)] = t(RMSEmat1[16:20,] < RMSEmat2[16:20,])
boldmat[, seq(2, 10, 2)] = t(RMSEmat1[16:20,] > RMSEmat2[16:20,])
boldmat[is.na(boldmat)] = FALSE
printbold(RMSExtable, which = boldmat, NA.string = "-")

############################################
###### PCT stfit better than gapfill #######
############################################
mat12 = RMSEmat1 < RMSEmat2
mat12
pctmat = matrix(NA, 4,3)
for(i in 1:4){
  for(j in 1:3){
    tmpvec = c(mat12[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    pctmat[i,j] = sum(tmpvec, na.rm=TRUE)/sum(!is.na(tmpvec))
  }
}
pctmat = t(pctmat)
rownames(pctmat) = c("< 10%", "40%~60%", "80%~90%")
colnames(pctmat) = c("S1", "S2", "S3", "S4")
xtable::xtable(pctmat)


##########################################################
###### Average RMSE by missing pattern and seasons #######
##########################################################
library(feather)
library(dplyr)
df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 2000] = NA

#### partial missing image indexes with different missing percentage
pidx0.1 = c(66, 75, 348, 573, 605)
pidx0.4_0.6 = c(74, 156, 273, 285, 326)
pidx0.8_0.95 = c(112, 184, 318, 448, 508)
pidx = c(pidx0.1, pidx0.4_0.6, pidx0.8_0.95)
n.na = apply(mat0[pidx,], 1, function(x) sum(is.na(x)))
n.na /961 * 100

#### integrated RMSE matrix for stfit
RMSEmat1.int = matrix(NA, 4,3)
weight = rep(n.na, each = 5)
for(i in 1:4){
  for(j in 1:3){
    tmpvec = c(RMSEmat1[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    nonna.idx = !-is.na(tmpvec)
    RMSEmat1.int[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
  }
}
RMSEmat1.int

#### integrated RMSE matrix for gapfill
RMSEmat2.int = matrix(NA, 4,3)
weight = rep(n.na, each = 5)
for(i in 1:4){
  for(j in 1:3){
    tmpvec = c(RMSEmat2[seq((i-1)*5+1, i*5), seq((j-1)*5+1, j*5)])
    nonna.idx = !-is.na(tmpvec)
    RMSEmat2.int[i,j] = sum(tmpvec[nonna.idx]*weight[nonna.idx])/sum(weight[nonna.idx])
  }
}
RMSEmat2.int

RMSEmat.int = matrix(NA, 3, 8)
RMSEmat.int[, c(1,3,5,7)] = t(RMSEmat1.int)
RMSEmat.int[, c(2,4,6,8)] = t(RMSEmat2.int)
RMSEmat.int
rownames(RMSEmat.int) = c("< 10%", "40%~60%", "80%~90%")
colnames(RMSEmat.int) = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4")
RMSEmat.int


