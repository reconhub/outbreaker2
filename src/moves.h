
Rcpp::List cpp_move_mu(Rcpp::List param, Rcpp::List data, Rcpp::List config,
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue);

Rcpp::List cpp_move_pi(Rcpp::List param, Rcpp::List data, Rcpp::List config,
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue);

Rcpp::List cpp_move_eps(Rcpp::List param, Rcpp::List data, Rcpp::List config,
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue);

Rcpp::List cpp_move_lambda(Rcpp::List param, Rcpp::List data, Rcpp::List config,
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue);

Rcpp::List cpp_move_poisson_scale(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                           Rcpp::RObject custom_ll = R_NilValue,
                           Rcpp::RObject custom_prior = R_NilValue);

Rcpp::List cpp_move_sigma(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                                  Rcpp::RObject custom_ll = R_NilValue,
                                  Rcpp::RObject custom_prior = R_NilValue);
  
Rcpp::List cpp_move_t_inf(Rcpp::List param, Rcpp::List data,
			  Rcpp::RObject list_custom_ll = R_NilValue);

Rcpp::List cpp_move_alpha(Rcpp::List param, Rcpp::List data,
			  Rcpp::RObject list_custom_ll = R_NilValue);

Rcpp::List cpp_move_swap_cases(Rcpp::List param, Rcpp::List data,
			  Rcpp::RObject list_custom_ll = R_NilValue);

Rcpp::List cpp_move_kappa(Rcpp::List param, Rcpp::List data,
			  Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue);

