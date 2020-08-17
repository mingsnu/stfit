// Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// nbr
IntegerVector nbr(int ii, int nRow, int nCol, int dRow, int dCol);
RcppExport SEXP _stfit_nbr(SEXP iiSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP dRowSEXP, SEXP dColSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< int >::type dRow(dRowSEXP);
    Rcpp::traits::input_parameter< int >::type dCol(dColSEXP);
    rcpp_result_gen = Rcpp::wrap(nbr(ii, nRow, nCol, dRow, dCol));
    return rcpp_result_gen;
END_RCPP
}
// lc_cov1_
double lc_cov1_(const NumericMatrix& X, const NumericMatrix& W, int ii, int jj, int nRow, int nCol, NumericVector& pidx);
RcppExport SEXP _stfit_lc_cov1_(SEXP XSEXP, SEXP WSEXP, SEXP iiSEXP, SEXP jjSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP pidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type pidx(pidxSEXP);
    rcpp_result_gen = Rcpp::wrap(lc_cov1_(X, W, ii, jj, nRow, nCol, pidx));
    return rcpp_result_gen;
END_RCPP
}
// sparse_emp_cov_est
DataFrame sparse_emp_cov_est(NumericMatrix X, int nRow, int nCol, int nnr);
RcppExport SEXP _stfit_sparse_emp_cov_est(SEXP XSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP nnrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< int >::type nnr(nnrSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_emp_cov_est(X, nRow, nCol, nnr));
    return rcpp_result_gen;
END_RCPP
}
// sparse_emp_cov_est1
DataFrame sparse_emp_cov_est1(NumericMatrix X, int nRow, int nCol, int nnr, NumericVector pidx);
RcppExport SEXP _stfit_sparse_emp_cov_est1(SEXP XSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP nnrSEXP, SEXP pidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< int >::type nnr(nnrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pidx(pidxSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_emp_cov_est1(X, nRow, nCol, nnr, pidx));
    return rcpp_result_gen;
END_RCPP
}
// sparse_lc_cov_est
DataFrame sparse_lc_cov_est(NumericMatrix X, NumericMatrix W, int nRow, int nCol, int nnr);
RcppExport SEXP _stfit_sparse_lc_cov_est(SEXP XSEXP, SEXP WSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP nnrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< int >::type nnr(nnrSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_lc_cov_est(X, W, nRow, nCol, nnr));
    return rcpp_result_gen;
END_RCPP
}
// sparse_lc_cov_est1
DataFrame sparse_lc_cov_est1(NumericMatrix X, NumericMatrix W, int nRow, int nCol, int nnr, NumericVector pidx);
RcppExport SEXP _stfit_sparse_lc_cov_est1(SEXP XSEXP, SEXP WSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP nnrSEXP, SEXP pidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< int >::type nnr(nnrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pidx(pidxSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_lc_cov_est1(X, W, nRow, nCol, nnr, pidx));
    return rcpp_result_gen;
END_RCPP
}
// lc_cov_1d
double lc_cov_1d(const NumericVector& ids, const NumericVector& time, const NumericVector& resid, const NumericVector& W, int t1, int t2);
RcppExport SEXP _stfit_lc_cov_1d(SEXP idsSEXP, SEXP timeSEXP, SEXP residSEXP, SEXP WSEXP, SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type resid(residSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(lc_cov_1d(ids, time, resid, W, t1, t2));
    return rcpp_result_gen;
END_RCPP
}
// lc_cov_1d_est
NumericMatrix lc_cov_1d_est(const NumericVector& ids, const NumericVector& time, const NumericVector& resid, const NumericVector& W, const NumericVector& tt);
RcppExport SEXP _stfit_lc_cov_1d_est(SEXP idsSEXP, SEXP timeSEXP, SEXP residSEXP, SEXP WSEXP, SEXP ttSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type resid(residSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type tt(ttSEXP);
    rcpp_result_gen = Rcpp::wrap(lc_cov_1d_est(ids, time, resid, W, tt));
    return rcpp_result_gen;
END_RCPP
}
// mean_est
NumericVector mean_est(NumericMatrix X, int nRow, int nCol, NumericMatrix W);
RcppExport SEXP _stfit_mean_est(SEXP XSEXP, SEXP nRowSEXP, SEXP nColSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nRow(nRowSEXP);
    Rcpp::traits::input_parameter< int >::type nCol(nColSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_est(X, nRow, nCol, W));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stfit_nbr", (DL_FUNC) &_stfit_nbr, 5},
    {"_stfit_lc_cov1_", (DL_FUNC) &_stfit_lc_cov1_, 7},
    {"_stfit_sparse_emp_cov_est", (DL_FUNC) &_stfit_sparse_emp_cov_est, 4},
    {"_stfit_sparse_emp_cov_est1", (DL_FUNC) &_stfit_sparse_emp_cov_est1, 5},
    {"_stfit_sparse_lc_cov_est", (DL_FUNC) &_stfit_sparse_lc_cov_est, 5},
    {"_stfit_sparse_lc_cov_est1", (DL_FUNC) &_stfit_sparse_lc_cov_est1, 6},
    {"_stfit_lc_cov_1d", (DL_FUNC) &_stfit_lc_cov_1d, 6},
    {"_stfit_lc_cov_1d_est", (DL_FUNC) &_stfit_lc_cov_1d_est, 5},
    {"_stfit_mean_est", (DL_FUNC) &_stfit_mean_est, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_stfit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
