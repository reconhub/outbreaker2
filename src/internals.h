
std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i);

size_t cpp_sample1(Rcpp::IntegerVector x);

size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i);

Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i);

Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i);

Rcpp::List cpp_swap_cases(Rcpp::List param, size_t i);

size_t cpp_get_n_mutations(Rcpp::List data, size_t i, size_t j);
