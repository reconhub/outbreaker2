
void copy_values(Rcpp::IntegerVector a, Rcpp::IntegerVector b);

size_t cpp_sample1(Rcpp::IntegerVector x);

size_t cpp_propose_alpha(Rcpp::IntegerVector t_inf, size_t i );

Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i);


