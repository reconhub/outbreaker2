#ifndef OUTBREAKER2_LIKELIHOODS_H
#define OUTBREAKER2_LIKELIHOODS_H

#include <Rcpp.h>

// Core likelihood functions

// [[Rcpp::export(rng = false)]]
double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
		      Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, size_t i,
		      Rcpp::RObject custom_function = R_NilValue);






// [[Rcpp::export(rng = false)]]
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
				Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, size_t i,
				Rcpp::RObject custom_function = R_NilValue);






// [[Rcpp::export(rng = false)]]
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
			      Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, size_t i,
			      Rcpp::RObject custom_function = R_NilValue);






// [[Rcpp::export(rng = false)]]
double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
			Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, size_t i,
			Rcpp::RObject custom_function = R_NilValue);






// [[Rcpp::export(rng = false)]]
double cpp_ll_contact(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
		      Rcpp::RObject custom_function = R_NilValue);

double cpp_ll_contact(Rcpp::List data, Rcpp::List param, size_t i,
		      Rcpp::RObject custom_function = R_NilValue);

// [[Rcpp::export(rng = false)]]
double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
		     Rcpp::RObject custom_functions = R_NilValue);

double cpp_ll_timing(Rcpp::List data, Rcpp::List param, size_t i,
		     Rcpp::RObject custom_functions = R_NilValue);

// [[Rcpp::export(rng = false)]]
double cpp_ll_patient_transfer(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
                               Rcpp::RObject custom_function = R_NilValue) 
  
double cpp_ll_patient_transfer(Rcpp::List data, Rcpp::List param, size_t i,
                                 Rcpp::RObject custom_function = R_NilValue)
  
  // [[Rcpp::export(rng = false)]]
double cpp_ll_potential_colonised(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
                                  Rcpp::RObject custom_function = R_NilValue) 

double cpp_ll_potential_colonised(Rcpp::List data, Rcpp::List param, size_t i,
                                    Rcpp::RObject custom_function = R_NilValue)

// Aggregated functions, i.e. summing some of the above


// [[Rcpp::export(rng = false)]]
double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i = R_NilValue,
		  Rcpp::RObject custom_functions = R_NilValue);

double cpp_ll_all(Rcpp::List data, Rcpp::List param, size_t i,
		  Rcpp::RObject custom_functions = R_NilValue);


#endif
