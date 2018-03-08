#' DOY mean estimation
#'
#' @param rst  a *Raster object
#' @param doy.idx a vector of DOY index of the *Raster object
#'
#' @return a *Raster object
#' @export
#'
#' @examples
doyMeanEst <- function(rst, doy.idx = 1:nlayers(rst)){
  .smoothFun <- function(x){
    nonna.idx = !is.na(x)
    splfit <- smooth.spline(doy.idx[nonna.idx], x[nonna.idx])
    predict(splfit, doy.idx)$y
  }
  calc(rst, .smoothFun)
}

#' Remove outlier
#'
#' An outlier is defined as points outside the whiskers of the boxplot
#' over the time domain (DOY).
#'
#' @param rst a *Raster object
#'
#' @return a *Raster object
#' @export
#'
#' @examples 
rmOutlier <- function(rst){
  .rmOutlier <- function(x){
    na.idx = is.na(x)
    boxstat = boxplot(x[!na.idx], plot=FALSE)
    x[!na.idx][(x[!na.idx] > boxstat$stats[5,]) | (x[!na.idx] < boxstat$stats[1,])] = NA
    x
  }
  calc(rst, .rmOutlier)
}

#' Image imputation using time of the day information
#'
#' @param rst.list A list of `RasterBrick` or `RasterStack` objects, whose elements are
#' measurements of the same location at different time of the day of the same day of 
#' the year
#'
#' @return
#' @export
#'
#' @examples
tod_impute <- function(rst.list, sim.mat){
  n.rst.list = length(rst.list)
  if(!ident(lapply(rst.list, dim)))
    stop("temporal_impute: elemetns of rst.list should have the same dimension")
  n.doy = nlayers(rst.list[[1]])
  doy.idx = 1:n.doy ## doy index start from 1
  mat.list = list()
  ## load data into matrix; each row is the obervation of a pixel over time (doy)
  for(i in 1:n.rst.list){
    mat.list[[i]] = matrix(as.integer(values(rst.list[[i]])), ncol = n.doy) ## save space using integer
  }
  if(missing(sim.mat))
    sim.mat = similarity_mat(mat.list)
  mat.imputed.list = mat.list
  ## format(object.size(mat.list), units = "Mb")
  ## if a pixel has no observation over the time, treat it as a "mask pixel" and remove it
  ## from imputation procedure
  mask.pixel.idx = which(apply(mat.imputed.list[[1]], 1, function(x) sum(!is.na(x))) == 0)
  ## loop through all pixel and do smoothing; smooth seperatedly for each time points of the day
  for(i in setdiff(1:nrow(mat.imputed.list[[1]]), mask.pixel.idx)){ ## the i-th pixel
    splfit.list = list()
    doy.missidx.mat = matrix(0, nrow = n.doy, ncol = n.rst.list)
    for(j in 1:n.rst.list){ ## the j-th time point (image) of the day
      doy.missidx.mat[,j] <- ifelse(is.na(mat.list[[j]][i, ]), 0, 1)
      nonna.idx <- which(doy.missidx.mat[,j] == 1) ## idx for non-na doy
      # splfit.list[[j]] <- smooth.spline(doy.idx[nonna.idx], mat.list[[j]][i, ][nonna.idx])
      # ## impute missing values using the smoothed curve (using the information of the sequence
      # ## of images at one single time point; some values will be updated later using information from
      # ## other observed time points of the same doy.
      # mat.imputed.list[[j]][i,][doy.idx[-nonna.idx]] = predict(splfit.list[[j]], doy.idx[-nonna.idx])$y
    }
    tmp = apply(doy.missidx.mat, 1, sum)
    fullobs.doy.idx = which(tmp == 4) ## fully observed doy index
    partmiss.doy.idx = which(tmp < 4)
    ## partmiss.doy.idx = which(tmp > 0 & tmp < 4) ## partial missing doy index
    ## allmiss.doy.idx = which(tmp == 0) ## all missing doy index
    # update imputed values for partially missing doy index images
    partmissidx.mat = doy.missidx.mat[partmiss.doy.idx,]
    
    for(j in 1:n.rst.list){
      j.miss.idx = which(partmissidx.mat[,j] == 0)
      j.partmissidx.mat = partmissidx.mat[j.miss.idx, ]
      for(k in sim.mat[,j]){ ## search the most similar time point one by one
        idx = which(j.partmissidx.mat[,k] == 1)
        tmp.doy.idx = partmiss.doy.idx[j.miss.idx[idx]]
        ## the esimated difference between j and k time points in tmp.doy.idx using the days
        ## before/after the current doy.
        n.tmp.doy.idx = length(tmp.doy.idx)
        tmp.doy.idx1 = tmp.doy.idx - 1
        tmp.doy.idx2 = tmp.doy.idx + 1
        if(tmp.doy.idx1[1] == 0) tmp.doy.idx1[1] = 2
        if(tmp.doy.idx2[n.tmp.doy.idx] == n.doy+1) tmp.doy.idx2[n.tmp.doy.idx] = n.doy - 1
        jk.diff = 0.5*(mat.imputed.list[[j]][i,][tmp.doy.idx1] - mat.imputed.list[[k]][i,][tmp.doy.idx1] +
                         mat.imputed.list[[j]][i,][tmp.doy.idx2] - mat.imputed.list[[k]][i,][tmp.doy.idx2])
        mat.imputed.list[[j]][i,][tmp.doy.idx] = mat.imputed.list[[k]][i,][tmp.doy.idx] + jk.diff
        ## after imputing with the k-th time point images, update the index
        j.miss.idx = j.miss.idx[-idx]
        j.partmissidx.mat = partmissidx.mat[j.miss.idx, ]
      }
    }
  }
  class(mat.imputed.list) = "GFtod"
  attr(mat.imputed.list, "rstinfo") = rst.list
  mat.imputed.list
}

#' Missing value percentages
#'
#' @param x A `RasterStack` object
#' @param mc.cores Numer of cores to use
#'
#' @return A vector of percent of missing values for each layer
#' @export
#'
#' @examples
pctMissing <- function(x, mc.cores){
  if(missing(mc.cores)) mc.cores = parallel::detectCores()
  doParallel::registerDoParallel(cores=mc.cores)
  nl = nlayers(x)
  foreach::foreach(i=1:nl, .combine = c) %dopar%{
    val = values(x[[i]])
    sum(is.na(val))/length(val)
  }
}

## test whether all input elemetns are equal
ident <- function(...) {
  args <- c(...)
  if (length(args) > 2L) {
    #  recursively call ident()
    out <- c(identical(args[1] , args[2]) , ident(args[-1]))
  } else{
    out <- identical(args[1] , args[2])
  }
  return(all(out))
}


#' get similarity matrix
#'
#' @param mat.list a list of matrix
#'
#' @return a matrix consists of the index of the matrices. For example, supose there are three matrices in the
#' list and the first column is c(3,2) which means matrix 3 is more similar to matrix 1 than matrix 2.
#' @export
#'
#' @examples
similarity_mat <- function(mat.list){
  mat = cbind(c(3, 2, 4), c(4, 1, 3), c(1, 2, 4), c(2, 1, 3))
}

#' Get missing layer index
#'
#' @param rst.list a RasterStack or RasterBrick object or a list of them
#'
#' @return index of the missing layers
#' @export
#'
#' @examples
getMissingLayers <- function(rst.list){
  if(inherits(rst.list, "RasterStackBrick"))
    return(which(is.infinite(rst.list@data@min) | is.na(rst.list@data@min)))
  if(is.list(rst.list))
    return(lapply(1:length(rst.list), function(i) 
      which(is.infinite(rst.list[[i]]@data@min) | is.na(rst.list[[i]]@data@min))))
}

#' Get the missing pattern (mask) of x
#' @name getMask
#' @rdname getMask
#' @exportMethod getMask
setGeneric (
  name = "getMask",
  def = function(object, ...) {
    standardGeneric("getMask")
  }
)

#' @rdname getMask
#' @aliases getMask,matrix-method
setMethod("getMask", "matrix",
  definition = function(object, tol = 0.95, ...) {
    return(
      apply(object, 2, function(x) {
      sum(is.na(x))/length(x) >= tol
    }))
  }
)
