#ifndef SCALE_H
#define SCALE_H

#include <Rcpp.h>
using namespace Rcpp;

double normal_consistency_constants(const std::size_t);
double const_c2(const std::size_t);
double median_cpp(NumericVector);
double Mscale(NumericVector, const double, const double);
double tau_scale(NumericVector, const double, const double);
NumericVector wni(NumericVector, const double, const double, const double);
NumericVector weight_factor(NumericVector, IntegerVector);
NumericMatrix new_centers(NumericMatrix, NumericVector, IntegerVector,
                          const std::size_t, NumericVector);

#endif