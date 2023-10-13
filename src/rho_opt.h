#ifndef RHO_OPT_H
#define RHO_OPT_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector rhoOpt(NumericVector, const double);
NumericVector psiOpt(NumericVector, const double);
NumericVector derpsiOpt(NumericVector, const double);

#endif