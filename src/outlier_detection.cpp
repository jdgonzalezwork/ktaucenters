#include "scale.h"
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//' Flag outliers
//'
//' @param cutoff quantile of chi-square to be used as a threshold for outliers
//' detection.
//' @param b break down point.
//' @param ktau ktaucenters results.
//'
//' @return
//' Numeric vector with the weight factor for each observation
//'
// [[Rcpp::export(".flag_outliers")]]
List flag_outliers(const double cutoff, const double b, List ktau) {
  NumericVector distances = ktau["di"];
  NumericVector squared_distances = pow(distances, 2);
  NumericMatrix centers = ktau["centers"];
  const std::size_t p = centers.cols();
  double thr = R::qchisq(cutoff, p, true, false);
  double c = normal_consistency_constants(p);
  double robust_scale = Mscale(distances, c, b);
  IntegerVector indices = seq_along(distances);
  ktau["outliers"] = indices[squared_distances > thr * pow(robust_scale, 2)];
  return ktau;
}