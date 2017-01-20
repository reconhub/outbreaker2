
Rcpp::List cpp_move_mu(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
		       Rcpp::RObject custom_prior);

Rcpp::List cpp_move_pi(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
		       Rcpp::RObject custom_prior);


Rcpp::List cpp_move_t_inf(Rcpp::List data, Rcpp::List param);

Rcpp::List cpp_move_alpha(Rcpp::List data, Rcpp::List param);

Rcpp::List cpp_move_swap_cases(Rcpp::List data, Rcpp::List param);

Rcpp::List cpp_move_kappa(Rcpp::List data, Rcpp::List param, 
			  Rcpp::List config);

