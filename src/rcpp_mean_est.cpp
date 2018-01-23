#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

double sumKernel(
    const NumericMatrix& X,    /* naip image */
  const NumericMatrix& W,    /* pre computed spatial weights */
  size_t i,      /* current location in rows */
  size_t j,      /* current location in columns */
  size_t nRow,   /* number of Rows */
  size_t nCol    /* number of Columns */
) {
  size_t dRow = W.nrow();
  size_t dCol = W.ncol();
  
  /* adjustment that must be applied for edge effects */
  size_t k, l, n;
  
  size_t k_start;
  size_t k_stop;
  size_t l_start;
  size_t l_stop;
  
  double sumYK = 0, sumK = 0;
  
  size_t k_local;
  size_t l_local;
  
  /* the starts */
  if( i < dRow/2 ) {
    k_start = 0; 
  } else {
    k_start = i - dRow/2 ;
  }
  if( j < dCol/2 ) {
    l_start = 0; 
  } else {
    l_start = j - dCol/2 ;
  }
  /* the stops */
  if( i + dRow/2 + 1 > nRow ) {
    k_stop = nRow; 
  } else {
    k_stop = i + dRow/2 + 1;
  }
  if( j + dCol/2 + 1  > nCol ) {
    l_stop = nCol; 
  } else {
    l_stop = j + dCol/2 + 1;
  }
  for(n = 0; n < X.nrow(); n++){
    for(k=k_start, k_local=k_start - i + (dRow/2); 
        k < k_stop; k++, k_local++) {
      for(l=l_start, l_local=l_start -j + (dCol/2);
          l < l_stop; l++, l_local++) {
        if(NumericVector::is_na(X(n, k * nCol + l))) continue;
        sumYK += X(n, k * nCol + l) * W(k_local, l_local);
        sumK += W(k_local, l_local);
      }
    }
  }
  if(sumK == 0.0){
    return NA_REAL;
  } else{
    return sumYK/sumK;
  }
}


// [[Rcpp::export]]
NumericVector mean_est(NumericMatrix X, size_t nRow, size_t nCol, NumericMatrix W) {
  // each row of X is a row stacked image
  // X.ncol() == nRow * nCol
  size_t i,j;
  NumericVector mu(nRow*nCol);
  
  for( i=0; i < nRow; i++) {
    for( j=0; j < nCol; j++) {
      mu[i*nCol + j] = sumKernel(X, W, i, j, nRow, nCol); 
    }
  }
  return mu;
}