#include <Rmath.h>
#include <Rcpp.h>
#include "internals.h"
#include "likelihoods.h"


// IMPORTANT: ON INDEXING VECTORS AND ANCESTRIES

// Most of the functions implemented here are susceptible to be called from R
// via Rcpp, and are therefore treated as interfaces. This causes a number of
// headaches when using indices of cases defined in R (1:N) to refer to elements
// in Rcpp / Cpp vectors (0:N-1). By convention, we store all data on the
// original scale (1:N), and modify indices whenever accessing elements of
// vectors. In other words, in an expression like 'alpha[j]', 'j' should always
// be on the internal scale (0:N-1).

// In all these functions, 'SEXP i' is an optional vector of case indices, on
// the 1:N scale.






// ---------------------------

// This likelihood corresponds to the probability of observing a number of
// mutations between cases and their ancestors. See src/likelihoods.cpp for
// details of the Rcpp implmentation.

// The likelihood is based on the number of mutations between a case and its
// ancestor; these are extracted from a pairwise genetic distance matrix
// (data$D) the log-likelihood is computed as: sum(mu^nmut + (1-mu)^(L-nmut))
// with:

// 'mu' is the mutation probability
// 'L' the number of sites in the alignment
// 'n_mut' the number of mutations between an ancestor and its descendent
// 'n_non_mut' the number of sites that have not mutated

// For any given case at 'nmut' mutations from its ancestor, with kappa
// generations in between, the log-likelihood is defined as:

// log(mu) * n_mut + log(1 - mu) * {(L - n_mut) + (L * (kappa-1))}


// when summing over several individuals, it becomes:

// log(mu) * sum_i(n_mut_i) + log(1-mu) * sum_i((L - n_mut_i) + (L * (kappa_i - 1)))

double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i,
		      Rcpp::RObject custom_function) {
  Rcpp::IntegerMatrix D = data["D"];
  if (D.ncol() < 1) return 0.0;

  size_t N = static_cast<size_t>(data["N"]);
  if (N < 2) return 0.0;

  if (custom_function == R_NilValue) {

    // Variables from the data & param
    Rcpp::NumericMatrix w_dens = data["log_w_dens"];
    size_t K = w_dens.nrow();
    double mu = Rcpp::as<double>(param["mu"]);
    long int L = Rcpp::as<int>(data["L"]);
    Rcpp::IntegerVector alpha = param["alpha"]; // values are on 1:N
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::LogicalVector has_dna = data["has_dna"];


    // Local variables used for computatoins
    size_t n_mut = 0;
    size_t n_non_mut = 0;
    double out = 0;
    bool found[1];
    size_t ances[1];
    size_t n_generations[1];
    found[0] = false;
    ances[0] = NA_INTEGER;
    n_generations[0] = NA_INTEGER;

  
    // Invalid values of mu
    if (mu < 0.0 || mu > 1.0) {
      return R_NegInf;
    }

    
    // NOTE ON MISSING SEQUENCES

    // Terms participating to the genetic likelihood correspond to pairs
    // of ancestor-descendent which have a genetic sequence. The
    // log-likelihood of other pairs is 0.0, and can therefore be
    // ommitted. Note the possible source of confusion in indices here:

    // 'has_dna' is a vector, thus indexed from 0:(N-1)
	  
    // 'cpp_get_n_mutations' is a function, and thus takes indices on 1:N


    
    // all cases are retained
    
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) { // 'j' on 0:(N-1)
	if (alpha[j] != NA_INTEGER) {

	  // kappa restriction

	  if (kappa[j] < 1 || kappa[j] > K) {
	    return R_NegInf;
	  }

	  // missing sequences handled here
	  
	  if (has_dna[j]) {
	 
	    lookup_sequenced_ancestor(alpha, kappa, has_dna, j + 1,
				      ances, n_generations, found);

	    if (found[0]) {

	      n_mut = cpp_get_n_mutations(data, j + 1, ances[0]); // remember the offset
	      n_non_mut = L - n_mut;

	      out += n_mut*log(n_generations[0]*mu) +
		n_non_mut*log(1 - n_generations[0]*mu);
		
	    }
	  }
	}
      }

    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
	size_t j = vec_i[k] - 1; // offset
	if (alpha[j] != NA_INTEGER) {
	  // kappa restriction
	  if (kappa[j] < 1 || kappa[j] > K) {
	    return R_NegInf;
	  }

	  // missing sequences handled here
	  	  
	  if (has_dna[j]) {
	 
	    lookup_sequenced_ancestor(alpha, kappa, has_dna, j + 1, 
				      ances, n_generations, found);

	    if (found[0]) {

	      n_mut = cpp_get_n_mutations(data, j + 1, ances[0]); // remember the offset
	      n_non_mut = L - n_mut;
	      
	      out += n_mut*log(n_generations[0]*mu) +
		n_non_mut*log(1 - n_generations[0]*mu);
	      
	    }
	  }
	  
	}

      }
    }

    return(out);
    
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(data, param));
  }
}


double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_genetic(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}





// ---------------------------

// This likelihood corresponds to the probability of observing infection dates
// of cases given the infection dates of their ancestors.

double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i,
				Rcpp::RObject custom_function) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  if (custom_function == R_NilValue) {

    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::NumericMatrix w_dens = data["log_w_dens"];
    size_t K = w_dens.nrow();

    double out = 0.0;

    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
	if (alpha[j] != NA_INTEGER) {
	  size_t delay = t_inf[j] - t_inf[alpha[j] - 1]; // offset
	  if (delay < 1 || delay > w_dens.ncol()) {
	    return  R_NegInf;
	  }
	  if (kappa[j] < 1 || kappa[j] > K) {
	    return  R_NegInf;
	  }

	  out += w_dens(kappa[j] - 1, delay - 1);
	}
      }
    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
	size_t j = vec_i[k] - 1; // offset
	if (alpha[j] != NA_INTEGER) {
	  size_t delay = t_inf[j] - t_inf[alpha[j] - 1]; // offset
	  if (delay < 1 || delay > w_dens.ncol()) {
	    return  R_NegInf;
	  }
	  if (kappa[j] < 1 || kappa[j] > K) {
	    return  R_NegInf;
	  }

	  out += w_dens(kappa[j] - 1, delay - 1);
	}

      }
    }

    return out;
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(data, param));
  }
}


double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timing_infections(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}





// ---------------------------

// This likelihood corresponds to the probability of reporting dates of cases
// given their infection dates.

double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i,
			      Rcpp::RObject custom_function) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  if (custom_function == R_NilValue) {

    Rcpp::IntegerVector dates = data["dates"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::NumericVector f_dens = data["log_f_dens"];

    double out = 0.0;

    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
	size_t delay = dates[j] - t_inf[j];
	if (delay < 1 || delay > f_dens.size()) {
	  return  R_NegInf;
	}
	out += f_dens[delay - 1];
      }
    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
	size_t j = vec_i[k] - 1; // offset
	size_t delay = dates[j] - t_inf[j];
	if (delay < 1 || delay > f_dens.size()) {
	  return  R_NegInf;
	}
	out += f_dens[delay - 1];
      }
    }

    return out;
  }  else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(data, param));
  }
}


double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timing_sampling(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}





// ---------------------------

// This likelihood corresponds to the probability of a given number of
// unreported cases on an ancestry.

// The likelihood is given by a geometric distribution with probability 'pi'
// to report a case

// - 'kappa' is the number of generation between two successive cases
// - 'kappa-1' is the number of unreported cases

double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i,
			Rcpp::RObject custom_function) {
  Rcpp::NumericMatrix w_dens = data["log_w_dens"];
  size_t K = w_dens.nrow();

  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  double pi = static_cast<double>(param["pi"]);
  Rcpp::IntegerVector kappa = param["kappa"];

  // p(pi < 0) = p(pi > 1) = 0
  if (pi < 0.0 || pi > 1.0) {
    return R_NegInf;
  }

  if (custom_function == R_NilValue) {

    double out = 0.0;

    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
	if (kappa[j] != NA_INTEGER) {
	  if (kappa[j] < 1 || kappa[j] > K) {
	    return  R_NegInf;
	  }
	  out += R::dgeom(kappa[j] - 1.0, pi, 1); // first arg must be cast to double
	}
      }
    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
	size_t j = vec_i[k] - 1; // offset
	if (kappa[j] != NA_INTEGER) {
	  if (kappa[j] < 1 || kappa[j] > K) {
	    return  R_NegInf;
	  }
	  out += R::dgeom(kappa[j] - 1.0, pi, 1); // first arg must be cast to double
	}
      }
    }

    return out;
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(data, param));
  }
}


double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_reporting(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}





// ---------------------------

// This likelihood corresponds to the probability of observing a a reported
// contact between cases and their ancestors. See
// src/likelihoods.cpp for details of the Rcpp implmentation.

// The likelihood is based on the contact status between a case and its
// ancestor; this is extracted from a pairwise contact matrix (data$C), the
// log-likelihood is computed as:
// true_pos*eps + false_pos*eps*xi +
// false_neg*(1- eps) + true_neg*(1 - eps*xi)
//
// with:
// 'eps' is the contact reporting coverage
// 'lambda' is the non-infectious contact rate
// 'true_pos' is the number of contacts between transmission pairs
// 'false_pos' is the number of contact between non-transmission pairs
// 'false_neg' is the number of transmission pairs without contact
// 'true_neg' is the number of non-transmission pairs without contact

double cpp_ll_contact(Rcpp::List data, Rcpp::List param, SEXP i,
		       Rcpp::RObject custom_function) {
  Rcpp::NumericMatrix contacts = data["contacts"];
  if (contacts.ncol() < 1) return 0.0;

  size_t C_combn = static_cast<size_t>(data["C_combn"]);
  size_t C_nrow = static_cast<size_t>(data["C_nrow"]);

  size_t N = static_cast<size_t>(data["N"]);
  if (N < 2) return 0.0;

  if (custom_function == R_NilValue) {

    double eps = Rcpp::as<double>(param["eps"]);
    double lambda = Rcpp::as<double>(param["lambda"]);
    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector kappa = param["kappa"];

    size_t true_pos = 0;
    size_t false_pos = 0;
    size_t false_neg = 0;
    size_t true_neg = 0;
    size_t imports = 0;
    size_t unobsv_case = 0;

    // p(eps < 0 || lambda < 0) = 0
    if (eps < 0.0 || lambda < 0.0) {
      return R_NegInf;
    }

    // all cases are retained (currently no support for i subsetting)
    for (size_t j = 0; j < N; j++) {
      if (alpha[j] == NA_INTEGER) {
	imports += 1;
      } else if (kappa[j] > 1) {
	unobsv_case += 1;
      } else {
	true_pos += contacts(j, alpha[j] - 1); // offset
      }
    }

    false_pos = C_nrow - true_pos;
    false_neg = N - imports - unobsv_case - true_pos;
    true_neg = C_combn - true_pos - false_pos - false_neg;

    // deal with special case when lambda == 0 and eps == 1, to avoid log(0)
    if(lambda == 0.0) {
      if(false_pos > 0) {
	return R_NegInf;
      } else {
	log(eps) * (double) true_pos +
	  log(1 - eps) * (double) false_neg +
	  log(1 - eps*lambda) * (double) true_neg;
      }
    } else if(eps == 1.0) {
      if(false_neg > 0) {
	return R_NegInf;
      } else {
	return log(eps) * (double) true_pos +
	  log(eps*lambda) * (double) false_pos +
	  log(1 - eps*lambda) * (double) true_neg;
      }
    } else {
      return log(eps) * (double) true_pos +
	log(eps*lambda) * (double) false_pos +
	log(1 - eps) * (double) false_neg +
	log(1 - eps*lambda) * (double) true_neg;
    }
  } else { //use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(data, param));
  }
}


double cpp_ll_contact(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_contact(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}





// ---------------------------

// This likelihood corresponds to the sums of the separate timing likelihoods,
// which include:

// - p(infection dates): see function cpp_ll_timing_infections
// - p(collection dates): see function cpp_ll_timing_sampling

double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i,
		     Rcpp::RObject custom_functions) {

  if (custom_functions == R_NilValue) {
    return cpp_ll_timing_infections(data, param, i) +
      cpp_ll_timing_sampling(data, param, i);
  } else { // use of a customized likelihood functions
    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);
    return cpp_ll_timing_infections(data, param, i, list_functions["timing_infections"]) +
      cpp_ll_timing_sampling(data, param, i, list_functions["timing_sampling"]);

  }
}


double cpp_ll_timing(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timing(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}




// ---------------------------

// This likelihood corresponds to the sums of the separate likelihoods, which
// include:

// - p(infection dates): see function cpp_ll_timing_infections
// - p(collection dates): see function cpp_ll_timing_sampling
// - p(genetic diversity): see function cpp_ll_genetic
// - p(missing cases): see function cpp_ll_reporting
// - p(contact): see function cpp_ll_contact 

double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i,
		  Rcpp::RObject custom_functions) {

  if (custom_functions == R_NilValue) {

    return cpp_ll_timing_infections(data, param, i) +
      cpp_ll_timing_sampling(data, param, i) +
      cpp_ll_genetic(data, param, i) +
      cpp_ll_reporting(data, param, i) +
      cpp_ll_contact(data, param, i);

  }  else { // use of a customized likelihood functions
    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);

    return cpp_ll_timing_infections(data, param, i, list_functions["timing_infections"]) +
      cpp_ll_timing_sampling(data, param, i, list_functions["timing_sampling"]) +
      cpp_ll_genetic(data, param, i, list_functions["genetic"]) +
      cpp_ll_reporting(data, param, i, list_functions["reporting"]) +
      cpp_ll_contact(data, param, i, list_functions["contact"]);

  }
}


double cpp_ll_all(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_all(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}
