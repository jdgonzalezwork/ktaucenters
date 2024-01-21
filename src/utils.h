#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

IntegerVector top_index(NumericVector, int, bool);
double median_cpp(NumericVector);
double max_tolerance(NumericMatrix, NumericMatrix);

#endif