#include "distance.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Distance Matrix Computation
//'
//' Computes and returns the distance matrix using euclidean distance
//' measure to compute the distances between the rows of a
//' data matrix.
//'
//' @param x a numeric matrix.
//'
//' @return
//' A numeric matrix with the distances between the rows of a matrix.
//'
// [[Rcpp::export(.distance)]]
NumericMatrix distance(NumericMatrix x) {
  const std::size_t n = x.nrow();
  double d;
  NumericMatrix out(no_init(n, n));

  for (std::size_t i = 0; i < n - 1; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      d = sqrt(sum(pow(x.row(i) - x.row(j), 2.0)));
      out(i, i) = 0.0;
      out(j, i) = d;
      out(i, j) = d;
    }
  }
  return out;
}
