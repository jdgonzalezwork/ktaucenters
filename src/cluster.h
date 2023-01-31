#ifndef CLUSTER_H
#define CLUSTER_H

#include <Rcpp.h>
using namespace Rcpp;

List cluster_location(NumericMatrix, NumericMatrix);
IntegerVector cluster_counter(IntegerVector, const std::size_t);

#endif