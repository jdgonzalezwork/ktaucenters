#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

List dist_to_kNN(NumericMatrix D, const std::size_t k) {
  //' Find the k Nearest Neighbors
  //'
  //' \strong{Ties}: If the kth and the (k+1)th nearest neighbor are tied, then
  //' the neighbor found first is returned and the other one is ignored.
  //'
  //' @param D a distance matrix, which contains the distances between the rows
  //' of a matrix.
  //' @param k number of neighbors to find.
  //'
  //' @return A list with the following components:
  //' \item{dist }{a matrix with distances. }
  //' \item{id }{a matrix with indices of the k nearest neighbors. }
  //' \item{k }{number neighbors used. }

  const std::size_t n = D.nrow();

  IntegerMatrix id(no_init(n, k));
  NumericMatrix d(no_init(n, k));
  IntegerVector neighbors_id(no_init(k));

  for (std::size_t i = 0; i < n; ++i) {
    NumericVector column_distance = D.column(i);
    column_distance[i] = R_PosInf;
    neighbors_id = top_index(column_distance, k, false);
    id(i, _) = neighbors_id;
    for (std::size_t j = 0; j < k; ++j) {
      d(i, j) = column_distance[neighbors_id[j]];
    }
  }
  return (List::create(_["dist"] = d, _["id"] = id, _["k"] = k));
}
