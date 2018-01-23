#include <Rcpp.h>
#include <stdlib.h>
#include <vector>
using namespace Rcpp;

std::vector<int> nbr_(int ii, int nRow, int nCol, int nnr){
  // (i, j) is the current point coordinates
  // bvec: to the Bottom of (i,j) vector
  // rvec: to the Right of (i,j) vector
  int k, l;
  std::vector<int> bvec, rvec, vec;
  int i = ii % nRow;
  int j = ii / nRow ;
  // Rcpp::Rcout << "i = " << (int)i << std::endl;
  // Rcpp::Rcout << "j = " << (int)j << std::endl;
  for(k = i; k <= std::min(i + nnr, nRow-1); k++){
    bvec.push_back(j*nRow + k);
  }
  
  if(j < nCol - 1){
    for(l = j+1; l <= std::min(j + nnr, nCol-1); l++){
      for(k = std::max(i - nnr, (int)0); k <= std::min(i + nnr, nRow-1); k++){
        rvec.push_back(l*nRow + k);
      }
    }
    vec.reserve( bvec.size() + rvec.size() ); // prebvecllocbvecte memory
    vec.insert( vec.end(), bvec.begin(), bvec.end() );
    vec.insert( vec.end(), rvec.begin(), rvec.end() );
    return vec;
  } else {
    return bvec;
  }
}

double emp_cov_(const NumericMatrix &X, int i, int j){
  double tmp, res = 0.0;
  int n = 0;
  for(int k = 0; k < X.nrow(); k++){
    tmp = X(k, i) * X(k, j);
    if(!NumericVector::is_na(tmp)){
      res += tmp;
      n += 1;
    }
  }
  if(n < 2)
    return NA_REAL;
  else
    return res/(n-1);
}

// [[Rcpp::export]]
DataFrame sparse_cov_est(NumericMatrix X, int nRow, int nCol, int nnr) {
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  std::vector<int> ridx;
  std::vector<int> cidx;
  std::vector<int> ii_nbr;
  std::vector<double> value;

  int ii, jj;
  for(ii = 0; ii < nRow*nCol; ii++){
    ii_nbr = nbr_(ii, nRow, nCol, nnr);
    for(jj = 0; jj < ii_nbr.size(); jj++){
      ridx.push_back(ii + 1);
      cidx.push_back(ii_nbr[jj] + 1);
      value.push_back(emp_cov_(X, ii, ii_nbr[jj]));
    }
  }
  return DataFrame::create( 
    _["ridx"]  = ridx, 
    _["cidx"]  = cidx, 
    _["value"] = value
  );
}
