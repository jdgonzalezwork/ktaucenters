#include "scale.h"
#include "cluster.h"
#include "rho_opt.h"
#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Normal Consistency Constants
//'
//' M scale tuning constants so it is consistent with the standard
//' normal distribution for the quasi optimal \eqn{\rho} function used in
//' \code{\link{rhoOpt}}. These constants were computed for 1 \eqn{\leq} p
//' \eqn{\leq} 400.
//'
//' @param p dimension where observation lives.
//'
//' @return tuning constant.
//'
//' @references
//' [1] Maronna, R. A., Martin, R. D., Yohai, V. J., & Salibián-Barrera, M.
//' (2018). 'Robust statistics: theory and methods (with ' R). Wiley.
//'
//' [2] Salibian-Barrera, M., Willems, G., & Zamar, R. (2008). The fast-tau
//' estimator for regression. 'Journal of Computational and Graphical
//' Statistics, 17(3), 659-682.
//'
//'@export
// [[Rcpp::export]]
double normal_consistency_constants(const std::size_t p) {
  // TODO: Replace hard coded values with function
  NumericVector vaux{
      0.404629, 0.6944748, 0.8985921, 1.063144, 1.204321, 1.329791, 1.443817,
      1.548994, 1.647149,  1.7395370, 1.827075, 1.910406, 1.990170, 2.066772,
      2.140529, 2.211772,  2.2807420, 2.347639, 2.412622, 2.475882, 2.537545,
      2.597723, 2.656494,  2.7140160, 2.770276, 2.825434, 2.879547, 2.932612,
      2.984741, 3.035955,  3.086320,  3.135869, 3.184648, 3.232684, 3.279986,
      3.326633, 3.372634,  3.418005,  3.462781, 3.506981, 3.550627, 3.593741,
      3.636342, 3.678449,  3.720075,  3.761236, 3.801950, 3.842231, 3.882080,
      3.921557, 3.960609,  3.999286,  4.037581, 4.075532, 4.113132, 4.150390,
      4.187306, 4.223932,  4.260211,  4.296178, 4.331845, 4.367245, 4.402347,
      4.437173, 4.471721,  4.506001,  4.540023, 4.573796, 4.607320, 4.640600,
      4.673648, 4.706454,  4.739037,  4.771396, 4.803549, 4.835481, 4.867181,
      4.898696, 4.930009,  4.961124,  4.992046, 5.022776, 5.053319, 5.083678,
      5.113838, 5.143877,  5.173678,  5.203352, 5.232814, 5.262137, 5.291296,
      5.320304, 5.349150,  5.377833,  5.406382, 5.434777, 5.463004, 5.491069,
      5.519042, 5.546857,  5.574513,  5.602029, 5.629455, 5.656713, 5.683843,
      5.710848, 5.737726,  5.764478,  5.791107, 5.817614, 5.844001, 5.870269,
      5.896420, 5.922455,  5.948373,  5.974178, 5.999901, 6.025461, 6.050950,
      6.076328, 6.101599,  6.126763,  6.151823, 6.176781, 6.201638, 6.226395,
      6.251055, 6.275624,  6.300084,  6.324457, 6.348727, 6.372907, 6.397029,
      6.421037, 6.444946,  6.468764,  6.492492, 6.516168, 6.539737, 6.563222,
      6.586623, 6.609943,  6.633180,  6.656336, 6.679371, 6.702367, 6.725285,
      6.748124, 6.770886,  6.793581,  6.816188, 6.838730, 6.861204, 6.883556,
      6.905896, 6.928129,  6.950294,  6.972394, 6.994428, 7.016405, 7.038284,
      7.060111, 7.081862,  7.103564,  7.125190, 7.146749, 7.168252, 7.189681,
      7.211047, 7.232349,  7.253587,  7.274746, 7.295884, 7.316926, 7.337911,
      7.358844, 7.379737,  7.400556,  7.421316, 7.442017, 7.462661, 7.483248,
      7.503778, 7.524252,  7.544670,  7.565039, 7.585339, 7.605592, 7.625797,
      7.645945, 7.666039,  7.686080,  7.706069, 7.726007, 7.745910, 7.765748,
      7.785535, 7.805273,  7.824960,  7.844558, 7.864147, 7.883689, 7.903179,
      7.922623, 7.942019,  7.961367,  7.980669, 7.999924, 8.019132, 8.038295,
      8.057414, 8.076486,  8.095513,  8.114495, 8.133433, 8.152326, 8.171176,
      8.189983, 8.208746,  8.227467,  8.246149, 8.264782, 8.283357, 8.301944,
      8.320429, 8.338918,  8.357330,  8.375742, 8.394071, 8.412383, 8.430632,
      8.448894, 8.467070,  8.485221,  8.503334, 8.521408, 8.539443, 8.557449,
      8.575411, 8.593294,  8.611180,  8.629047, 8.646855, 8.664627, 8.682362,
      8.700059, 8.717721,  8.735375,  8.752975, 8.770500, 8.788031, 8.805526,
      8.822969, 8.840417,  8.857788,  8.875119, 8.892458, 8.909744, 8.926998,
      8.944220, 8.961409,  8.978539,  8.995658, 9.012743, 9.029797, 9.046818,
      9.063807, 9.080764,  9.097691,  9.114586, 9.131449, 9.148281, 9.165083,
      9.181853, 9.198593,  9.215303,  9.231982, 9.248631, 9.265251, 9.281840,
      9.298400, 9.314931,  9.331432,  9.347904, 9.364348, 9.380762, 9.397148,
      9.413505, 9.429834,  9.446134,  9.462407, 9.478651, 9.494868, 9.511057,
      9.527219, 9.543353,  9.559460,  9.575540, 9.591593, 9.607619, 9.623618,
      9.639591, 9.655538,  9.671458,  9.687352, 9.703220, 9.719062, 9.734879,
      9.750670, 9.766435,  9.782175,  9.797895, 9.813593, 9.829224, 9.844870,
      9.860490, 9.876084,  9.891652,  9.907196, 9.922714, 9.938208, 9.953677,
      9.969121, 9.984542,  9.999938,  10.01531, 10.03066, 10.04598, 10.06129,
      10.07656, 10.09182,  10.10705,  10.12226, 10.13745, 10.15261, 10.16776,
      10.18288, 10.19798,  10.21304,  10.22809, 10.24312, 10.25813, 10.27312,
      10.28810, 10.30303,  10.31794,  10.33285, 10.34773, 10.36260, 10.37743,
      10.39222, 10.40702,  10.42181,  10.43656, 10.45129, 10.46599, 10.48068,
      10.49537, 10.51000,  10.52463,  10.53923, 10.55382, 10.56839, 10.58294,
      10.59746, 10.61197,  10.62646,  10.64092, 10.65537, 10.66980, 10.68421,
      10.69860, 10.71297,  10.72732,  10.74165, 10.75597, 10.77026, 10.78454,
      10.79879, 10.81303,  10.82725,  10.84145, 10.85564, 10.86980, 10.88395,
      10.89807, 10.91218,  10.92627,  10.94030, 10.95436, 10.96841, 10.98243,
      10.99643, 11.01041,  11.02437,  11.03832, 11.05225, 11.06617, 11.08006,
      11.09393};
  return vaux(p - 1);
}

double const_c2(const std::size_t p) {
  //' Tuning constant for \eqn{\tau}-estimator
  //'
  //' Tuning constant approximation for \eqn{\tau}-estimator to reach 90%
  //' efficiency
  //'
  //' @param p dimension where observation lives.
  //'
  //' @return
  //' Tuning constant for \eqn{\tau}-estimator.
  //'
  //' @references
  //' Maronna, R. A. and Yohai, V. J. (2017). Robust and efficient estimation of
  //' multivariate scatter and location. Computational Statistics & Data
  //' Analysis, 109:64–75.

  return 2.9987 * pow(p, -0.4647);
}

//' M scale
//'
//' The M scale of an univariate sample.
//'
//' @param u numeric vector with positive values.
//' @param c a tuning constant. If consistency to standard normal distribution
//' is desired use \code{\link{normal_consistency_constants}}.
//' @param b the desired break down point.
//'
//' @return
//' M scale value.
//'
//' @examples
//' Mscale(u = rnorm(100), c = 1, b = 0.5)
//'
//' @references Maronna, R. A., Martin, R. D., Yohai, V. J., & Salibian-Barrera,
//' M. (2018). Robust statistics: theory and methods (with R). Wiley.
//'
//'@export
// [[Rcpp::export]]
double Mscale(NumericVector u, const double c, const double b) {

  const double max_error_diff = 1e-10;

  // TODO: Indicate why dividing by 0.6745
  double sn = median_cpp(abs(u)) / 0.6745;
  if (sn == 0)
    return sn;

  double diff = mean(rhoOpt(u / sn, c)) - b;
  if (fabs(diff) <= max_error_diff)
    return sn;

  while (diff > 0.0) {
    sn = 1.5 * sn;
    diff = mean(rhoOpt(u / sn, c)) - b;
  }

  unsigned int i = 0;
  double error = 1.0;
  double Fk, Gk, factor;
  NumericVector var(no_init(u.size()));
  while ((i < 1000) & (error > max_error_diff)) {
    var = u / sn;
    Fk = mean(rhoOpt(var, c));
    Gk = mean(psiOpt(var, c) * var);
    factor = (Fk - Gk - b) / (2.0 * Fk - Gk - 2.0 * b);
    error = fabs(factor - 1);
    sn = sn * fabs(factor);
    i += 1;
  }
  return sn;
}

double tau_scale(NumericVector u, const double c, const double s) {
  //' \eqn{\tau} scale
  //'
  //' The \eqn{\tau} scale of an univariate sample.
  //'
  //' @param u numeric vector with positive values.
  //' @param s M scale for the same univariate sample (u).
  //' @param c tuning constant to regulate estimator's efficiency.
  //'
  //' @return
  //' \eqn{\tau} scale value

  return s * sqrt(mean(rhoOpt(u / s, c)));
}

NumericVector wni(NumericVector distances, const double c1, const double c2,
                  const double s) {

  //' Weight function
  //'
  //' @param distances numeric vector with distances from each point to its
  //' cluster center.
  //' @param c1 tuning constant used for m scale estimation.
  //' @param c2 tuning constant used for \eqn{\tau} scale estimation.
  //' @param s M scale for the same vector of distances.
  //'
  //' @return
  //' Numeric vector with the weight for each observation

  NumericVector dnor = distances / s;
  const double A = sum(2 * rhoOpt(dnor, c2) - psiOpt(dnor, c2) * dnor);
  const double B = sum(psiOpt(dnor, c1) * dnor);
  return ifelse(dnor == 0.0,
                (A / (3.25 * pow(c1, 2))) + (B / (3.25 * pow(c2, 2))),
                (A * psiOpt(dnor, c1) + B * psiOpt(dnor, c2)) / dnor);
}

NumericVector weight_factor(NumericVector wni, IntegerVector clusters) {

  //' Weight factor
  //'
  //' @param wni numeric vector with the weight for each observation.
  //' @param clusters integer vector with the cluster location for each
  //' observation.
  //'
  //' @return
  //' Numeric vector with the weight factor for each observation

  if (wni.size() != clusters.size()) {
    stop("Both wni and clusters must have same size.");
  }

  const std::size_t n_records = wni.size();
  const std::size_t n_cluster = unique(clusters).size();

  NumericVector out(no_init(n_records));
  NumericVector wni_sum(n_cluster);
  IntegerVector cluster_count(n_cluster);

  for (std::size_t idx = 0; idx < n_records; ++idx) {
    wni_sum[clusters[idx] - 1] += wni[idx];
    cluster_count[clusters[idx] - 1] += 1;
  }

  for (std::size_t idx = 0; idx < n_records; ++idx) {
    if (wni_sum[clusters[idx] - 1] != 0.0) {
      out[idx] = wni[idx] / wni_sum[clusters[idx] - 1];
    } else {
      out[idx] = 1.0 / cluster_count[clusters[idx] - 1];
    }
  }

  return out;
}

NumericMatrix new_centers(NumericMatrix x, NumericVector weights,
                          IntegerVector cluster_location,
                          const std::size_t n_clusters,
                          NumericVector distances) {

  //' Computes new cluster centers
  //'
  //' @param x numeric matrix of size n x p with all observations.
  //' @param weights numeric vector of size n with the weight factor for each
  //' observation.
  //' @param cluster_location integer vector of size n with the cluster location
  //' for each observation.
  //' @param n_clusters integer vector of size n with the
  //' cluster location for each ' observation.
  //' @param distances numeric vector of size n with the distances from each
  //' observation to its nearest cluster center.
  //'
  //' @return
  //' Numeric matrix with the new cluster centers.

  const std::size_t p = x.cols();

  NumericMatrix out(n_clusters, p);

  for (std::size_t column = 0; column != p; ++column) {
    NumericVector weighted_column = x.column(column) * weights;
    NumericVector sums(n_clusters, 0.0);

    for (int idx = 0; idx != cluster_location.size(); ++idx) {
      sums[cluster_location[idx] - 1] += weighted_column[idx];
    }

    for (auto &cluster_it : cluster_location) {
      out(cluster_it - 1, column) = sums[cluster_it - 1];
    }
  }

  LogicalVector empty_clusters =
      cluster_counter(cluster_location, n_clusters) == 0;

  if (any(empty_clusters).is_true()) {
    std::vector<int> empty_pos;
    for (int i = 0; i < empty_clusters.size(); ++i) {
      if (empty_clusters[i]) {
        empty_pos.push_back(i);
      }
    }

    IntegerVector furthest_indices =
        top_index(distances, empty_pos.size(), true);

    for (std::size_t column = 0; column != p; ++column) {
      std::size_t idx = 0;
      for (auto &position : empty_pos) {
        out(position, column) = x(furthest_indices[idx], column);
        ++idx;
      }
    }
  }

  return out;
}
