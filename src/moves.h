
Rcpp::List cpp_move_mu(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
		       Rcpp::RObject custom_ll = R_NilValue, 
		       Rcpp::RObject custom_prior = R_NilValue);

Rcpp::List cpp_move_pi(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
		       Rcpp::RObject custom_ll = R_NilValue, 
		       Rcpp::RObject custom_prior = R_NilValue);


Rcpp::List cpp_move_t_inf(Rcpp::List data, Rcpp::List param,
			  Rcpp::RObject list_custom_ll = R_NilValue);

Rcpp::List cpp_move_alpha(Rcpp::List data, Rcpp::List param,
			  Rcpp::RObject list_custom_ll = R_NilValue);

Rcpp::List cpp_move_swap_cases(Rcpp::List data, Rcpp::List param,
			  Rcpp::RObject list_custom_ll = R_NilValue);

Rcpp::List cpp_move_kappa(Rcpp::List data, Rcpp::List param, 
			  Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue);

