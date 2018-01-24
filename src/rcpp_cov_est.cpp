#include <Rcpp.h>
#include <stdlib.h>
#include <vector>
using namespace Rcpp;

std::vector<size_t> nbr_(size_t ii, size_t nRow, size_t nCol, size_t nnr){
    // The "following" neighorhood pixel indexes (including ii)
    // (i, j) is the current posize_t coordinates
    // bvec: to the Bottom of (i,j) vector
    // rvec: to the Right of (i,j) vector
    size_t k, l; //cursor coordinates
    std::vector<size_t> bvec, rvec, vec;
    size_t i = ii / nCol;
    size_t j = ii % nCol;
    // Rcpp::Rcout << "i = " << (size_t)i << std::endl;
    // Rcpp::Rcout << "j = " << (size_t)j << std::endl;
    for(l = j; l <= std::min(j + nnr, nCol-1); l++){
      rvec.push_back(i*nCol + l);
    }
    
    if(i < nRow - 1){
      for(k = i+1; k <= std::min(i + nnr, nRow-1); k++){
        for(l = std::max(j - nnr, (size_t)0); l <= std::min(j + nnr, nCol-1); l++){
          bvec.push_back(k*nCol + l);
        }
      }
      vec.reserve( rvec.size() + bvec.size()); // prebvecllocbvecte memory
      vec.insert( vec.end(), rvec.begin(), rvec.end() );
      vec.insert( vec.end(), bvec.begin(), bvec.end() );
      return vec;
    } else {
      return rvec;
    }
}


std::vector<size_t> nbr_(size_t ii, size_t nRow, size_t nCol, size_t dRow, size_t dCol){
  // All neiborhood pixel indexes (including ii) within a given dRow x dCol rectangle
  // where ii is the centroid of the rectangle.
  // nRow, nCol are the data matrix dimension
  // dRow, dCol are the Weight matrix dimension.
  size_t k, l;
  size_t i = ii / nCol;
  size_t j = ii % nCol;
  std::vector<size_t> vec;
  for(k = std::max(i - dRow/2, (size_t)0); k <= std::min(i + dRow/2, nRow-1); k++){
    for(l = std::max(j - dCol/2, (size_t)0); l <= std::min(j + dCol/2, nCol-1); l++){
      vec.push_back(k*nCol + l);
    }
  }
  return vec;
}

double emp_cov_(const NumericMatrix &X, size_t i, size_t j){
  //sparse empirical covariance estimation
  double tmp, res = 0.0;
  size_t n = 0;
  for(size_t k = 0; k < X.nrow(); k++){
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
DataFrame sparse_cov_est(NumericMatrix X, size_t nRow, size_t nCol, size_t nnr) {
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  std::vector<size_t> ridx;
  std::vector<size_t> cidx;
  std::vector<size_t> ii_nbr;
  std::vector<double> value;

  size_t ii, jj;
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
