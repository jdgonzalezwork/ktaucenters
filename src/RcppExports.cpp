// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cluster_location
List cluster_location(NumericMatrix x, NumericMatrix centers);
RcppExport SEXP _ktaucenters_cluster_location(SEXP xSEXP, SEXP centersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    rcpp_result_gen = Rcpp::wrap(cluster_location(x, centers));
    return rcpp_result_gen;
END_RCPP
}
// cluster_counter
IntegerVector cluster_counter(IntegerVector cluster_location, const std::size_t n_clusters);
RcppExport SEXP _ktaucenters_cluster_counter(SEXP cluster_locationSEXP, SEXP n_clustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type cluster_location(cluster_locationSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type n_clusters(n_clustersSEXP);
    rcpp_result_gen = Rcpp::wrap(cluster_counter(cluster_location, n_clusters));
    return rcpp_result_gen;
END_RCPP
}
// distance
NumericMatrix distance(NumericMatrix x);
RcppExport SEXP _ktaucenters_distance(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(distance(x));
    return rcpp_result_gen;
END_RCPP
}
// dist_to_kNN
List dist_to_kNN(NumericMatrix D, const std::size_t k);
RcppExport SEXP _ktaucenters_dist_to_kNN(SEXP DSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_to_kNN(D, k));
    return rcpp_result_gen;
END_RCPP
}
// ktaucenters_run
List ktaucenters_run(NumericMatrix x, NumericMatrix centers, const double tolerance, const std::size_t max_iter);
RcppExport SEXP _ktaucenters_ktaucenters_run(SEXP xSEXP, SEXP centersSEXP, SEXP toleranceSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ktaucenters_run(x, centers, tolerance, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// flag_outliers
List flag_outliers(const double cutoff, const double b, List ktau);
RcppExport SEXP _ktaucenters_flag_outliers(SEXP cutoffSEXP, SEXP bSEXP, SEXP ktauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< List >::type ktau(ktauSEXP);
    rcpp_result_gen = Rcpp::wrap(flag_outliers(cutoff, b, ktau));
    return rcpp_result_gen;
END_RCPP
}
// rho_opt
NumericVector rho_opt(NumericVector x, const double c);
RcppExport SEXP _ktaucenters_rho_opt(SEXP xSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_opt(x, c));
    return rcpp_result_gen;
END_RCPP
}
// psi_opt
NumericVector psi_opt(NumericVector x, const double c);
RcppExport SEXP _ktaucenters_psi_opt(SEXP xSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_opt(x, c));
    return rcpp_result_gen;
END_RCPP
}
// derpsi_opt
NumericVector derpsi_opt(NumericVector x, double c);
RcppExport SEXP _ktaucenters_derpsi_opt(SEXP xSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(derpsi_opt(x, c));
    return rcpp_result_gen;
END_RCPP
}
// point_density
NumericVector point_density(NumericMatrix D, const std::size_t k);
RcppExport SEXP _ktaucenters_point_density(SEXP DSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(point_density(D, k));
    return rcpp_result_gen;
END_RCPP
}
// robin_center
std::size_t robin_center(NumericVector idp, IntegerVector indices, const double crit_robin);
RcppExport SEXP _ktaucenters_robin_center(SEXP idpSEXP, SEXP indicesSEXP, SEXP crit_robinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type idp(idpSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const double >::type crit_robin(crit_robinSEXP);
    rcpp_result_gen = Rcpp::wrap(robin_center(idp, indices, crit_robin));
    return rcpp_result_gen;
END_RCPP
}
// robinden
List robinden(NumericMatrix D, const std::size_t n_clusters, const std::size_t mp);
RcppExport SEXP _ktaucenters_robinden(SEXP DSEXP, SEXP n_clustersSEXP, SEXP mpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type n_clusters(n_clustersSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type mp(mpSEXP);
    rcpp_result_gen = Rcpp::wrap(robinden(D, n_clusters, mp));
    return rcpp_result_gen;
END_RCPP
}
// normal_consistency_constants
double normal_consistency_constants(const std::size_t p);
RcppExport SEXP _ktaucenters_normal_consistency_constants(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::size_t >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_consistency_constants(p));
    return rcpp_result_gen;
END_RCPP
}
// const_c2
double const_c2(const std::size_t p);
RcppExport SEXP _ktaucenters_const_c2(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::size_t >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(const_c2(p));
    return rcpp_result_gen;
END_RCPP
}
// mscale
double mscale(NumericVector u, const double c, const double b);
RcppExport SEXP _ktaucenters_mscale(SEXP uSEXP, SEXP cSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mscale(u, c, b));
    return rcpp_result_gen;
END_RCPP
}
// tau_scale
double tau_scale(NumericVector u, const double c, const double s);
RcppExport SEXP _ktaucenters_tau_scale(SEXP uSEXP, SEXP cSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(tau_scale(u, c, s));
    return rcpp_result_gen;
END_RCPP
}
// wni
NumericVector wni(NumericVector distances, const double c1, const double c2, const double s);
RcppExport SEXP _ktaucenters_wni(SEXP distancesSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< const double >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const double >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(wni(distances, c1, c2, s));
    return rcpp_result_gen;
END_RCPP
}
// weight_factor
NumericVector weight_factor(NumericVector wni, IntegerVector clusters);
RcppExport SEXP _ktaucenters_weight_factor(SEXP wniSEXP, SEXP clustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type wni(wniSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type clusters(clustersSEXP);
    rcpp_result_gen = Rcpp::wrap(weight_factor(wni, clusters));
    return rcpp_result_gen;
END_RCPP
}
// new_centers
NumericMatrix new_centers(NumericMatrix x, NumericVector weights, IntegerVector cluster_location, const std::size_t n_clusters, NumericVector distances);
RcppExport SEXP _ktaucenters_new_centers(SEXP xSEXP, SEXP weightsSEXP, SEXP cluster_locationSEXP, SEXP n_clustersSEXP, SEXP distancesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cluster_location(cluster_locationSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type n_clusters(n_clustersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distances(distancesSEXP);
    rcpp_result_gen = Rcpp::wrap(new_centers(x, weights, cluster_location, n_clusters, distances));
    return rcpp_result_gen;
END_RCPP
}
// median_cpp
double median_cpp(NumericVector x);
RcppExport SEXP _ktaucenters_median_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(median_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// top_index
IntegerVector top_index(NumericVector v, int n, bool decreasing);
RcppExport SEXP _ktaucenters_top_index(SEXP vSEXP, SEXP nSEXP, SEXP decreasingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type decreasing(decreasingSEXP);
    rcpp_result_gen = Rcpp::wrap(top_index(v, n, decreasing));
    return rcpp_result_gen;
END_RCPP
}
// max_tolerance
double max_tolerance(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _ktaucenters_max_tolerance(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(max_tolerance(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ktaucenters_cluster_location", (DL_FUNC) &_ktaucenters_cluster_location, 2},
    {"_ktaucenters_cluster_counter", (DL_FUNC) &_ktaucenters_cluster_counter, 2},
    {"_ktaucenters_distance", (DL_FUNC) &_ktaucenters_distance, 1},
    {"_ktaucenters_dist_to_kNN", (DL_FUNC) &_ktaucenters_dist_to_kNN, 2},
    {"_ktaucenters_ktaucenters_run", (DL_FUNC) &_ktaucenters_ktaucenters_run, 4},
    {"_ktaucenters_flag_outliers", (DL_FUNC) &_ktaucenters_flag_outliers, 3},
    {"_ktaucenters_rho_opt", (DL_FUNC) &_ktaucenters_rho_opt, 2},
    {"_ktaucenters_psi_opt", (DL_FUNC) &_ktaucenters_psi_opt, 2},
    {"_ktaucenters_derpsi_opt", (DL_FUNC) &_ktaucenters_derpsi_opt, 2},
    {"_ktaucenters_point_density", (DL_FUNC) &_ktaucenters_point_density, 2},
    {"_ktaucenters_robin_center", (DL_FUNC) &_ktaucenters_robin_center, 3},
    {"_ktaucenters_robinden", (DL_FUNC) &_ktaucenters_robinden, 3},
    {"_ktaucenters_normal_consistency_constants", (DL_FUNC) &_ktaucenters_normal_consistency_constants, 1},
    {"_ktaucenters_const_c2", (DL_FUNC) &_ktaucenters_const_c2, 1},
    {"_ktaucenters_mscale", (DL_FUNC) &_ktaucenters_mscale, 3},
    {"_ktaucenters_tau_scale", (DL_FUNC) &_ktaucenters_tau_scale, 3},
    {"_ktaucenters_wni", (DL_FUNC) &_ktaucenters_wni, 4},
    {"_ktaucenters_weight_factor", (DL_FUNC) &_ktaucenters_weight_factor, 2},
    {"_ktaucenters_new_centers", (DL_FUNC) &_ktaucenters_new_centers, 5},
    {"_ktaucenters_median_cpp", (DL_FUNC) &_ktaucenters_median_cpp, 1},
    {"_ktaucenters_top_index", (DL_FUNC) &_ktaucenters_top_index, 3},
    {"_ktaucenters_max_tolerance", (DL_FUNC) &_ktaucenters_max_tolerance, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ktaucenters(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
