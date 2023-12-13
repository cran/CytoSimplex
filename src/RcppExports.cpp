// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// colNormalize_dense
arma::mat colNormalize_dense(arma::mat x, arma::vec colsums);
RcppExport SEXP _CytoSimplex_colNormalize_dense(SEXP xSEXP, SEXP colsumsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type colsums(colsumsSEXP);
    rcpp_result_gen = Rcpp::wrap(colNormalize_dense(x, colsums));
    return rcpp_result_gen;
END_RCPP
}
// is_rawCounts_sparse
bool is_rawCounts_sparse(const arma::sp_mat& x);
RcppExport SEXP _CytoSimplex_is_rawCounts_sparse(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(is_rawCounts_sparse(x));
    return rcpp_result_gen;
END_RCPP
}
// is_rawCounts_dense
bool is_rawCounts_dense(const arma::mat& x);
RcppExport SEXP _CytoSimplex_is_rawCounts_dense(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(is_rawCounts_dense(x));
    return rcpp_result_gen;
END_RCPP
}
// euclidean_dense
arma::mat euclidean_dense(arma::mat& query, arma::mat& target);
RcppExport SEXP _CytoSimplex_euclidean_dense(SEXP querySEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type query(querySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(euclidean_dense(query, target));
    return rcpp_result_gen;
END_RCPP
}
// euclidean_sparse
arma::mat euclidean_sparse(arma::sp_mat query, arma::sp_mat target);
RcppExport SEXP _CytoSimplex_euclidean_sparse(SEXP querySEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type query(querySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(euclidean_sparse(query, target));
    return rcpp_result_gen;
END_RCPP
}
// cosine_dense
arma::mat cosine_dense(arma::mat query, arma::mat target);
RcppExport SEXP _CytoSimplex_cosine_dense(SEXP querySEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type query(querySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(cosine_dense(query, target));
    return rcpp_result_gen;
END_RCPP
}
// cosine_sparse
arma::mat cosine_sparse(arma::sp_mat query, arma::sp_mat target);
RcppExport SEXP _CytoSimplex_cosine_sparse(SEXP querySEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type query(querySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(cosine_sparse(query, target));
    return rcpp_result_gen;
END_RCPP
}
// rankMatrix_dense
Rcpp::List rankMatrix_dense(arma::mat& X);
RcppExport SEXP _CytoSimplex_rankMatrix_dense(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(rankMatrix_dense(X));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rank_matrix_dgc
std::vector<std::list<float> > cpp_rank_matrix_dgc(arma::vec& x, const arma::vec& p, int nrow, int ncol);
RcppExport SEXP _CytoSimplex_cpp_rank_matrix_dgc(SEXP xSEXP, SEXP pSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rank_matrix_dgc(x, p, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// rowAggregateSum_dense
arma::mat rowAggregateSum_dense(const arma::mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_rowAggregateSum_dense(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(rowAggregateSum_dense(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}
// rowAggregateSum_sparse
arma::mat rowAggregateSum_sparse(arma::sp_mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_rowAggregateSum_sparse(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(rowAggregateSum_sparse(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}
// colAggregateSum_dense
arma::mat colAggregateSum_dense(const arma::mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_colAggregateSum_dense(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(colAggregateSum_dense(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}
// colAggregateSum_sparse
arma::mat colAggregateSum_sparse(arma::sp_mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_colAggregateSum_sparse(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(colAggregateSum_sparse(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}
// colNNZAggr_dense
arma::mat colNNZAggr_dense(const arma::mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_colNNZAggr_dense(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(colNNZAggr_dense(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}
// colNNZAggr_sparse
arma::mat colNNZAggr_sparse(arma::sp_mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_colNNZAggr_sparse(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(colNNZAggr_sparse(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}
// rowNNZAggr_sparse
arma::mat rowNNZAggr_sparse(arma::sp_mat& X, const arma::uvec& groups, unsigned ngroups);
RcppExport SEXP _CytoSimplex_rowNNZAggr_sparse(SEXP XSEXP, SEXP groupsSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(rowNNZAggr_sparse(X, groups, ngroups));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CytoSimplex_colNormalize_dense", (DL_FUNC) &_CytoSimplex_colNormalize_dense, 2},
    {"_CytoSimplex_is_rawCounts_sparse", (DL_FUNC) &_CytoSimplex_is_rawCounts_sparse, 1},
    {"_CytoSimplex_is_rawCounts_dense", (DL_FUNC) &_CytoSimplex_is_rawCounts_dense, 1},
    {"_CytoSimplex_euclidean_dense", (DL_FUNC) &_CytoSimplex_euclidean_dense, 2},
    {"_CytoSimplex_euclidean_sparse", (DL_FUNC) &_CytoSimplex_euclidean_sparse, 2},
    {"_CytoSimplex_cosine_dense", (DL_FUNC) &_CytoSimplex_cosine_dense, 2},
    {"_CytoSimplex_cosine_sparse", (DL_FUNC) &_CytoSimplex_cosine_sparse, 2},
    {"_CytoSimplex_rankMatrix_dense", (DL_FUNC) &_CytoSimplex_rankMatrix_dense, 1},
    {"_CytoSimplex_cpp_rank_matrix_dgc", (DL_FUNC) &_CytoSimplex_cpp_rank_matrix_dgc, 4},
    {"_CytoSimplex_rowAggregateSum_dense", (DL_FUNC) &_CytoSimplex_rowAggregateSum_dense, 3},
    {"_CytoSimplex_rowAggregateSum_sparse", (DL_FUNC) &_CytoSimplex_rowAggregateSum_sparse, 3},
    {"_CytoSimplex_colAggregateSum_dense", (DL_FUNC) &_CytoSimplex_colAggregateSum_dense, 3},
    {"_CytoSimplex_colAggregateSum_sparse", (DL_FUNC) &_CytoSimplex_colAggregateSum_sparse, 3},
    {"_CytoSimplex_colNNZAggr_dense", (DL_FUNC) &_CytoSimplex_colNNZAggr_dense, 3},
    {"_CytoSimplex_colNNZAggr_sparse", (DL_FUNC) &_CytoSimplex_colNNZAggr_sparse, 3},
    {"_CytoSimplex_rowNNZAggr_sparse", (DL_FUNC) &_CytoSimplex_rowNNZAggr_sparse, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CytoSimplex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}