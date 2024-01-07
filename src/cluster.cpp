#include "cluster.h"
#include <Rcpp.h>
using namespace Rcpp;

List cluster_location(NumericMatrix x, NumericMatrix centers) {

  //' Returns each observation's cluster location based on the euclidean
  // distance
  //' between the observation and cluster centers.
  //'
  //' @param x a numeric matrix with observations.
  //' @param centers a numeric matrix with observations with cluster's centers
  //' coordinates.
  //'
  //' @return A list with the following components:
  //' \item{clusters }{a vector indicating each point observation's cluster. }
  //' \item{distance }{a vector with the distance from each
  //' observation to its corresponding cluster center. }

  const std::size_t k = centers.rows();
  const std::size_t n = x.rows();
  const std::size_t p = x.cols();
  IntegerVector clusters(no_init(n));
  NumericVector distance(no_init(n));
  double best, dd, tmp;
  std::size_t inew = 0;

  for (std::size_t n_iter = 0; n_iter < n; ++n_iter) {
    best = R_PosInf;
    for (std::size_t k_iter = 0; k_iter < k; ++k_iter) {
      dd = 0.0;
      for (std::size_t p_iter = 0; p_iter < p; ++p_iter) {
        tmp = x[n_iter + n * p_iter] - centers[k_iter + k * p_iter];
        dd += tmp * tmp;
      }
      if (dd < best) {
        best = dd;
        inew = k_iter + 1; // Base 1 cluster location.
      }
    }
    distance[n_iter] = std::sqrt(best);
    clusters[n_iter] = inew;
  }
  return List::create(_["clusters"] = clusters, _["distance"] = distance);
}

IntegerVector cluster_counter(IntegerVector cluster_location,
                              const int n_clusters) {

  //' Counts the number of observations per cluster
  //'
  //' @param cluster_location a vector with positive intergers containing each
  //' point observation's cluster.
  //' @param n_clusters total number of clusters
  //'
  //' @return
  //' Vector with number of observations per cluster

  IntegerVector counter(n_clusters);
  for (const auto &obs_cluster : cluster_location) {
    if (obs_cluster > 0 && obs_cluster <= n_clusters)
      counter[obs_cluster - 1] += 1;
  }
  return counter;
}
