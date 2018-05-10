library(xtable)
RMSEmat1 = readRDS("./pidx0.1/RMSEmat1.rds")
RMSEmat2 = readRDS("./pidx0.1/RMSEmat2.rds")
tmp1 = readRDS("./pidx0.4_0.6/RMSEmat1.rds")
tmp2 = readRDS("./pidx0.4_0.6/RMSEmat2.rds")
RMSEmat1 = cbind(RMSEmat1, tmp1)
RMSEmat2 = cbind(RMSEmat2, tmp2)
tmp1 = readRDS("./pidx0.8_0.95/RMSEmat1.rds")
tmp2 = readRDS("./pidx0.8_0.95/RMSEmat2.rds")
RMSEmat1 = cbind(RMSEmat1, tmp1)
RMSEmat2 = cbind(RMSEmat2, tmp2)

rownames(RMSEmat1) = paste0("F", 1:6)
colnames(RMSEmat1) = paste0("P", 1:15)
xtable::xtable(RMSEmat1)

rownames(RMSEmat2) = paste0("F", 1:6)
colnames(RMSEmat2) = paste0("P", 1:15)
xtable::xtable(RMSEmat2)

## combine
RMSEmat = matrix(0, 15, 12)
RMSEmat[, seq(1,12, 2)] = t(RMSEmat1)
RMSEmat[, seq(2,12, 2)] = t(RMSEmat2)
RMSEmat
rownames(RMSEmat) = paste0("P", 1:15)
colnames(RMSEmat) = paste0("F", rep(1:6, each=2))


RMSExtable = xtable::xtable(RMSEmat)

boldmat = matrix(TRUE, 15, 12)
## table(c(RMSEmat1 < RMSEmat2))
boldmat[, seq(1, 12, 2)] = t(RMSEmat1 < RMSEmat2)
boldmat[, seq(2, 12, 2)] = t(RMSEmat1 > RMSEmat2)
printbold(RMSExtable, which = boldmat)

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
