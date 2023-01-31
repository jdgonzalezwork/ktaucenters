#ifndef RHO_OPT_H
#define RHO_OPT_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector rho_opt(NumericVector, const double);
NumericVector psi_opt(NumericVector, const double);
NumericVector derpsi_opt(NumericVector, const double);

#endif