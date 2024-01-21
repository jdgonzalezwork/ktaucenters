#ifndef KNN_H
#define KNN_H

#include <Rcpp.h>
using namespace Rcpp;

List dist_to_kNN(NumericMatrix, const std::size_t);

#endif