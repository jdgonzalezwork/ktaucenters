#include "rho_opt.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Quasi optimal \eqn{\rho} function
//'
//' @param x numeric vector with positive values.
//' @param cc tunning constant.
//'
//' @return
//' Numeric vector with quasi optimal \eqn{\rho} computation for each element
//' of x.
//'
//'@references
//' [1] Salibian-Barrera, M., Willems, G., & Zamar, R. (2008). The fast-tau
//' estimator for regression. Journal of Computational and GraphicalStatistics,
//' 17(3), 659-682.
//'
//'@export
// [[Rcpp::export]]
NumericVector rhoOpt(NumericVector x, const double cc) {
  NumericVector out(no_init(x.size()));

  NumericVector::iterator out_it = out.begin();
  for (const auto &x_it : x) {
    if (fabs(x_it) <= 2 * cc) {
      *out_it = 0.5 * pow(x_it / cc, 2) / 3.25;
    } else if (fabs(x_it) <= 3 * cc) {
      *out_it = (1.792 - 0.972 * pow(x_it / cc, 2) + 0.432 * pow(x_it / cc, 4) -
                 0.052 * pow(x_it / cc, 6) + 0.002 * pow(x_it / cc, 8)) /
                3.25;
    } else {
      *out_it = 1.0;
    }
    ++out_it;
  }

  return out;
}

//' Derivative of the quasi optimal \eqn{\rho} function
//'
//' @param x numeric vector with positive values.
//' @param cc tunning constant.
//'
//' @return
//' Numeric vector with the derivative of the quasi optimal \eqn{\rho}
//' computation for each element of x.
//'
//'@export
// [[Rcpp::export]]
NumericVector psiOpt(NumericVector x, const double cc) {
  NumericVector out(no_init(x.size()));

  NumericVector::iterator out_it = out.begin();
  for (const auto &x_it : x) {
    if (fabs(x_it) <= 2 * cc) {
      *out_it = x_it / (3.25 * pow(cc, 2));
    } else if (fabs(x_it) <= 3 * cc) {
      *out_it =
          (-1.944 * x_it / pow(cc, 2) + 1.728 * pow(x_it, 3) / pow(cc, 4) -
           0.312 * pow(x_it, 5) / pow(cc, 6) +
           0.016 * pow(x_it, 7) / pow(cc, 8)) /
          3.25;
    } else {
      *out_it = 0.0;
    }
    ++out_it;
  }

  return out;
}

//' Second derivative of the quasi \eqn{\rho} function
//'
//' @param x numeric vector with positive values.
//' @param cc tunning constant.
//'
//' @return
//' Numeric vector with the second derivative of the quasi optimal \eqn{\rho}
//' computation for each element of x.
//'
//'@export
// [[Rcpp::export]]
NumericVector derpsiOpt(NumericVector x, double cc) {
  NumericVector out(no_init(x.size()));

  NumericVector::iterator out_it = out.begin();
  for (const auto &x_it : x) {
    if (fabs(x_it) <= 2 * cc) {
      *out_it = 1.0 / (3.25 * pow(cc, 2));
    } else if (fabs(x_it) <= 3 * cc) {
      *out_it = (-1.944 / pow(cc, 2) + 5.184 * pow(x_it, 2) / pow(cc, 4) -
                 1.56 * pow(x_it, 4) / pow(cc, 6) +
                 0.112 * pow(x_it, 6) / pow(cc, 8)) /
                3.25;
    } else {
      *out_it = 0.0;
    }
    ++out_it;
  }

  return out;
}
