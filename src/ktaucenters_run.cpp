#include "cluster.h"
#include "scale.h"
#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Robust Clustering algorithm based on centers, a robust and efficient version
//' of kmeans.
//'
//' @param x numeric matrix of size n x p with all observations.
//' @param centers numeric matrix with initial cluster centers.
//' @param tolerance maximum difference between current and new computed
//' clusters. Parameter used for the algorithm stopping rule.
//' @param max_iter a maximum number of iterations used for the algorithm
//' stopping rule.
//'
//' @return A list with the following components:
//' \item{tau }{\eqn{\tau} scale value. }
//' \item{iter }{number of iterations until convergence is achieved or maximum
//' number of iteration is reached. }
//' \item{di }{distance of each observation to its nearest cluster center. }
//' \item{centers }{numeric matrix of size K x p, with the estimated K centers.
//'  }
//' \item{clusters }{integer vector of size n with the cluster location for each
//' observation. }
//'
//'@references
//' [1] Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
//' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198.
//'
//' [2] Maronna, R. A. and Yohai, V. J. (2017). Robust and efficient estimation
//' of multivariate scatter and location.Computational Statistics &Data
//' Analysis, 109 : 64–75.
//'
// [[Rcpp::export(".ktaucenters_run")]]
List ktaucenters_run(NumericMatrix x, NumericMatrix centers,
                     const double tolerance, const std::size_t max_iter) {

  const std::size_t n_clusters = centers.rows();
  const std::size_t n = x.rows();
  const std::size_t p = x.cols();
  // According to Maronna and Yohai's reference
  const double c1 = 1.0;
  const double c2 = const_c2(p);
  // Target breakdown point
  const double b1 = 0.5;

  std::size_t iter = 0;
  double tol = tolerance + 1.0;
  NumericVector weights(n);
  NumericVector distance_min(n);
  NumericVector Wni(no_init(n));
  IntegerVector clusters(n);
  double tau, s;
  List cluster_loc;
  NumericMatrix old_centers(no_init(n_clusters, p));

  while (iter < max_iter && tol > tolerance) {
    // Step 1: (re)compute labels
    cluster_loc = cluster_location(x, centers);
    distance_min = cluster_loc["distance"];
    clusters = cluster_loc["clusters"];

    s = Mscale(distance_min, c1, b1);
    tau = tau_scale(distance_min, c2, s);
    // Step 2: (re)compute centers
    old_centers = centers;
    Wni = wni(distance_min, c1, c2, s);
    weights = weight_factor(Wni, clusters);
    centers = new_centers(x, weights, clusters, n_clusters, distance_min);

    tol = max_tolerance(old_centers, centers);
    iter += 1;
  }

  if (iter >= max_iter) {
    warning("Algorithm did not converge after %i iterations", max_iter);
  }

  return (List::create(_["centers"] = centers, _["cluster"] = clusters,
                       _["tau"] = tau, _["di"] = distance_min,
                       _["iter"] = iter));
}