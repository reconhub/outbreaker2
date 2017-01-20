

// Core likelihood functions

double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i,
		      Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, size_t i,
		      Rcpp::RObject custom_function = R_NilValue);






double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i,
				Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, size_t i,
				Rcpp::RObject custom_function = R_NilValue);






double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i,
			      Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, size_t i,
			      Rcpp::RObject custom_function = R_NilValue);






double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i,
			Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, size_t i,
			Rcpp::RObject custom_function = R_NilValue);






// Aggregated functions, i.e. summing some of the above

double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i,
		     Rcpp::RObject custom_functions = R_NilValue);

double cpp_ll_timing(Rcpp::List data, Rcpp::List param, size_t i,
		     Rcpp::RObject custom_functions = R_NilValue);






double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i,
		  Rcpp::RObject custom_functions = R_NilValue);

double cpp_ll_all(Rcpp::List data, Rcpp::List param, size_t i,
		  Rcpp::RObject custom_function = R_NilValue);
