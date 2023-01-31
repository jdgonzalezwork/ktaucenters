#include "cluster.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Returns each observation's cluster location based on the euclidean distance
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
//'
// [[Rcpp::export]]
List cluster_location(NumericMatrix x, NumericMatrix centers) {

  const std::size_t k = centers.rows();
  const std::size_t n = x.rows();
  IntegerVector clusters(no_init(n));
  NumericVector distance(no_init(n));
  NumericVector dist(k);

  for (std::size_t n_iter = 0; n_iter < n; ++n_iter) {
    for (std::size_t k_iter = 0; k_iter < k; ++k_iter) {
      dist[k_iter] = sqrt(sum(pow(x.row(n_iter) - centers.row(k_iter), 2.0)));
    }
    distance(n_iter) = min(dist);
    clusters(n_iter) = which_min(dist) + 1; // Base 1 cluster location.
  }
  return List::create(_["clusters"] = clusters, _["distance"] = distance);
}

//' Counts the number of observations per cluster
//'
//' @param cluster_location a vector with positive intergers containing each
//' point observation's cluster.
//' @param n_clusters total number of clusters
//'
//' @return
//' Vector with number of observations per cluster
//'
// [[Rcpp::export]]
IntegerVector cluster_counter(IntegerVector cluster_location,
                              const std::size_t n_clusters) {
  IntegerVector counter(n_clusters);
  for (const auto &obs_cluster : cluster_location) {
    if (obs_cluster > 0 && obs_cluster <= n_clusters)
      counter[obs_cluster - 1] += 1;
  }
  return counter;
}
