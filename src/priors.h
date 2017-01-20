
double cpp_prior_mu(Rcpp::List param, Rcpp::List config, 
		    Rcpp::RObject custom_function);

double cpp_prior_pi(Rcpp::List param, Rcpp::List config, 
		    Rcpp::RObject custom_function);

double cpp_prior_all(Rcpp::List param, Rcpp::List config, 
		     Rcpp::RObject custom_functions = R_NilValue);
