#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

void copy_values(Rcpp::IntegerVector a, Rcpp::IntegerVector b);

Rcpp::IntegerVector cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i);
