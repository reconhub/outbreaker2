#include <Rmath.h>
#include <Rcpp.h>
#include <algorithm>
#include "internals.h"
#include "likelihoods.h"


// IMPORTANT: ON INDEXING VECTORS AND ANCESTRIES

// Most of the functions implemented here are susceptible to be ed from R
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

  // Distance matrix
  Rcpp::IntegerMatrix D = data["D"];
  if (D.ncol() < 1) return 0.0;

  // Number of cases
  size_t N = static_cast<size_t>(data["N"]);
  if (N < 2) return 0.0;

  // Genetic sequence model
  std::string genetic_model = Rcpp::as<std::string>(data["genetic_model"]);

  Rcpp::List l;
  if (custom_function != R_NilValue) {
    l = Rcpp::as<Rcpp::List>(custom_function);
  }
  if (custom_function == R_NilValue || (l.size() > 0 && l[0] == R_NilValue)) {

    // Variables from the data & param
    Rcpp::NumericMatrix w_dens = data["log_w_dens"];
    size_t K = w_dens.nrow();
    double mu = Rcpp::as<double>(param["mu"]);
    long int L = Rcpp::as<int>(data["L"]);
    bool has_ctd_timed = data["has_ctd_timed"];
    Rcpp::IntegerVector alpha = param["alpha"]; // values are on 1:N
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector t_onw = param["t_onw"];
    Rcpp::IntegerVector t_sam = data["dna_dates"];
    Rcpp::IntegerVector id_in_dna = data["id_in_dna"];
    Rcpp::LogicalVector has_dna = data["has_dna"];

    // Invalid values of mu
    if (mu < 0.0 || mu > 1.0) {
      return R_NegInf;
    }

    // Initialise likelihood    
    double out = 0;
    
    if(genetic_model == "default") {

      // Local variables used for computations
      size_t n_mut = 0;
      size_t n_non_mut = 0;
      bool found[1];
      size_t ances[1];
      size_t n_generations[1];
      found[0] = false;
      ances[0] = NA_INTEGER;
      n_generations[0] = NA_INTEGER;
    
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

    } else if(genetic_model == "mrca") {

      Rcpp::NumericMatrix mrca = param["mrca"];
      Rcpp::NumericMatrix combn = data["dna_combn"];

      // add at the beginning of vector - indexing is now by case ID! (not ID - 1)
      t_inf.insert(t_inf.begin(), 100000);
      t_onw.insert(t_onw.begin(), 100000);

      int t_inf_1, t_inf_2;
      int id_1, id_2;
      int nmut, t_div, t_diff;
    
      //    double mut_sum = 0, time_sum = 0;

      // for(size_t j = 0; j < N; j++) {
      //   if (alpha[j] != NA_INTEGER) {
      for(size_t j = 0; j < mrca.nrow(); j++) {

	// get infection times of two nodes
	// we index by case ID because we have added a value to the beginning of t_inf
	id_1 = mrca(j, 0);
	id_2 = mrca(j, 1);

	// if we are imputing t_onw, use this infection time
	if(has_ctd_timed) {

	  // if we have an unobserved case, use the infection time from t_onw
	  if(kappa[id_1-1] == 1) {
	    t_inf_1 = t_inf[id_1];
	  } else {
	    t_inf_1 = t_onw[id_1];
	  }

	  if(kappa[id_2-1] == 1) {
	    t_inf_2 = t_inf[id_2];
	  } else {
	    t_inf_2 = t_onw[id_2];
	  }
	
	} else {
	
	  t_inf_1 = t_inf[id_1];
	  t_inf_2 = t_inf[id_2];
	  
	}

	// the earlier infection time represents the time of divergence. if one
	// case is the infector of the other (i.e. mrca = 0), then we choose the
	// other infection date (because we've set it to 10000)
	if(t_inf_1 > t_inf_2) {
	  t_div = t_inf_2;
	} else {
	  t_div = t_inf_1;
	}

	// get total divergence time (divergence to sampling)
	t_diff = abs(t_sam[id_in_dna[combn(j, 0)-1]-1] - t_div) +
	  abs(t_sam[id_in_dna[combn(j, 1)-1]-1] - t_div);

	// get number of mutations
	nmut = combn(j,2);

	//      printf("t_diff = %i | n_mut = %i \n", t_diff, nmut);
      
	// calculate likelihood
	out += R::dpois(nmut, t_diff*mu*L, true);
 
     }

    }

    return(out);
    
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(l[0]);
    int arity = l[1];
    if (arity == 3) return Rcpp::as<double>(f(data, param, i));
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

// This likelihood corresponds to the probability of observing
// infection dates of cases given the infection dates of their
// ancestors. It also imposes the strict assumption that *serial
// intervals* (i.e. delay between onset times) must be greater than
// zero.

double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i,
				Rcpp::RObject custom_function) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  Rcpp::List l;
  if (custom_function != R_NilValue) {
    l = Rcpp::as<Rcpp::List>(custom_function);
  }
  if (custom_function == R_NilValue || (l.size() > 0 && l[0] == R_NilValue)) {

    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::IntegerVector dates = data["dates"];
    Rcpp::NumericMatrix w_dens = data["log_w_dens"];

    bool has_ctd_timed = data["has_ctd_timed"];
    
    size_t K = w_dens.nrow();
    bool negative_si = data.containsElementNamed("negative_si") ? 
      Rcpp::as<bool>(data["negative_si"]) : true;

    double out = 0.0;
    double ll_1 = 0.0;
    double ll_2 = 0.0;

    size_t delay;
    size_t j;
    
    // Use the default temporal likelihood
    if(!has_ctd_timed) {
      // all cases are retained
      if (i == R_NilValue) {
	for (size_t j = 0; j < N; j++) {
	  if (alpha[j] != NA_INTEGER) {
	    // delay in infection times
	    int delay_inf = t_inf[j] - t_inf[alpha[j] - 1]; // offset
	    // delay in onset times
	    int delay_ons = dates[j] - dates[alpha[j] - 1]; // offset
	    if (delay_inf < 1 || delay_inf > w_dens.ncol()) {
	      return  R_NegInf;
	    }
	    if (!negative_si && delay_ons < 0) {
	      return  R_NegInf;
	    }
	    if (kappa[j] < 1 || kappa[j] > K) {
	      return  R_NegInf;
	    }
	    out += w_dens(kappa[j] - 1, delay_inf - 1);
	  }
	}
      } else {
	// only the cases listed in 'i' are retained
	size_t length_i = static_cast<size_t>(LENGTH(i));
	Rcpp::IntegerVector vec_i(i);
	for (size_t k = 0; k < length_i; k++) {
	  size_t j = vec_i[k] - 1; // offset
	  if (alpha[j] != NA_INTEGER) {
	    // delay in infection times
	    int delay_inf = t_inf[j] - t_inf[alpha[j] - 1]; // offset
	    // delay in onset times
	    int delay_ons = dates[j] - dates[alpha[j] - 1]; // offset
	    if (delay_inf < 1 || delay_inf > w_dens.ncol()) {
	      return  R_NegInf;
	    }
	    if (!negative_si && delay_ons < 0) {
	      return  R_NegInf;
	    }
	    if (kappa[j] < 1 || kappa[j] > K) {
	      return  R_NegInf;
	    }

	    out += w_dens(kappa[j] - 1, delay_inf - 1);
	  }
	}
      }
      // Infer the onward infection time
    } else {
      // all cases are retained
      Rcpp::IntegerVector t_onw = param["t_onw"];
      Rcpp::NumericMatrix w_unobs = data["log_w_unobs"];
      if (i == R_NilValue) {
	for (size_t j = 0; j < N; j++) {
	  if (alpha[j] != NA_INTEGER) {
	    if (kappa[j] < 1 || kappa[j] - 1 > K) {
	      return  R_NegInf;
	    }
	    // If we have no unobserved cases, use the default likelihood
	    if (kappa[j] == 1) {
	      delay = t_inf[j] - t_inf[alpha[j] - 1]; // offset
	      if (delay < 1 || delay > w_dens.ncol()) {
		return  R_NegInf;
	      }
	      out += w_dens(kappa[j] - 1, delay - 1);
	    } else if(kappa[j] > 1) {
	      // First account for the time between the infectee and the
	      // earliest unobserved case in the unobserved transmission chain -
	      // this uses w_unobs generation time (e.g. if the asymptomatic
	      // generation time is much longer)
	      delay = t_inf[j] - t_onw[j];
	      if (delay < 1 || delay > w_unobs.ncol()) {
		return  R_NegInf;
	      }

	      // We subtract two (one for indexing, one because we incorporate
	      // the second generation time in the next step
	      ll_1 = w_unobs(kappa[j] - 2, delay - 1);

	      // Account for time between infection of infector and earliest
	      // unobserved case - this follows the normal w_dens
	      delay = t_onw[j] - t_inf[alpha[j] - 1]; // offset
	      if (delay < 1 || delay > w_dens.ncol()) {
		return  R_NegInf;
	      }
	      // This is always going to be one generation (t_unobs is defined
	      // as such)

	      ll_2 = w_dens(0, delay - 1);

	      out += (ll_1 + ll_2)/2;
	      
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
	    if (kappa[j] < 1 || kappa[j] > K) {
	      return  R_NegInf;
	    }
	    // If we have no unobserved cases, use the default likelihood
	    if (kappa[j] == 1) {
	      size_t delay = t_inf[j] - t_inf[alpha[j] - 1]; // offset
	      if (delay < 1 || delay > w_dens.ncol()) {
		return  R_NegInf;
	      }
	      out += w_dens(kappa[j] - 1, delay - 1);
	    } else if(kappa[j] > 1) {
	      // First account for the time between the infectee and the
	      // earliest unobserved case in the unobserved transmission chain -
	      // this uses w_unobs generation time (e.g. if the asymptomatic
	      // generation time is much longer)
	      size_t delay = t_inf[j] - t_onw[j];
	      if (delay < 1 || delay > w_unobs.ncol()) {
		return  R_NegInf;
	      }
	      // We subtract two (one for indexing, one because we incorporate
	      // the second generation time in the next step
	      ll_1 = w_unobs(kappa[j] - 2, delay - 1);
	      // Account for time between infection of infector and earliest
	      // unobserved case - this follows the normal w_dens
	      delay = t_onw[j] - t_inf[alpha[j] - 1]; // offset
	      if (delay < 1 || delay > w_dens.ncol()) {
		return  R_NegInf;
	      }
	      // This is always going to be one generation (t_unobs is defined
	      // as such)
	      ll_2 = w_dens(0, delay - 1);

	      out += (ll_1 + ll_2)/2;
	      
	    }
	  }
	}
      }
    }

    return out;
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(l[0]);
    int arity = l[1];
    if (arity == 3) return Rcpp::as<double>(f(data, param, i));
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

  Rcpp::List l;
  if (custom_function != R_NilValue) {
    l = Rcpp::as<Rcpp::List>(custom_function);
  }
  if (custom_function == R_NilValue || (l.size() > 0 && l[0] == R_NilValue)) {

    Rcpp::IntegerVector dates = data["dates"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector id_in_f = data["id_in_f"];
    Rcpp::NumericMatrix f_dens = data["log_f_dens"];

    double out = 0.0;
    size_t delay;
    size_t j;

    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
	delay = dates[j] - t_inf[j];
	if (delay < 1 || delay > f_dens.nrow()) {
	  return  R_NegInf;
	}
	out += f_dens(delay - 1, id_in_f[j] - 1);
      }
    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
	j = vec_i[k] - 1; // offset
	delay = dates[j] - t_inf[j];
	if (delay < 1 || delay > f_dens.nrow()) {
	  return  R_NegInf;
	}
	out += f_dens(delay - 1, id_in_f[j] - 1);
      }
    }

    return out;
  }  else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(l[0]);
    int arity = l[1];
    if (arity == 3) return Rcpp::as<double>(f(data, param, i));
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

  Rcpp::List l;
  if (custom_function != R_NilValue) {
    l = Rcpp::as<Rcpp::List>(custom_function);
  }
  if (custom_function == R_NilValue || (l.size() > 0 && l[0] == R_NilValue)) {

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
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(l[0]);
    int arity = l[1];
    if (arity == 3) return Rcpp::as<double>(f(data, param, i));
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
// 'eta' is the contact sensitivity
// 'lambda' is the non-infectious contact rate
// 'true_pos' is the number of contacts between transmission pairs
// 'false_pos' is the number of contact between non-transmission pairs
// 'false_neg' is the number of transmission pairs without contact
// 'true_neg' is the number of non-transmission pairs without contact

double cpp_ll_contact(Rcpp::List data, Rcpp::List param, SEXP i,
		      Rcpp::RObject custom_function) {

  size_t N = static_cast<size_t>(data["N"]);
  if (N < 2) return 0.0;

  Rcpp::List l;
  if (custom_function != R_NilValue) {
    l = Rcpp::as<Rcpp::List>(custom_function);
  }
  if (custom_function == R_NilValue || (l.size() > 0 && l[0] == R_NilValue)) {

    Rcpp::List ctd_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_matrix"]);
    
    if(ctd_matrix_list.size() < 1) return 0.0;
    
    int C_ind;
    double out = 0;

    size_t C_combn = static_cast<size_t>(data["C_combn"]);
    Rcpp::IntegerVector C_nrow = Rcpp::as<Rcpp::IntegerVector>(data["C_nrow"]);

    Rcpp::NumericVector eps = param["eps"];
    Rcpp::NumericVector eta = param["eta"];
    Rcpp::NumericVector lambda = param["lambda"];
    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector t_onw = param["t_onw"];
    
    size_t list_size = ctd_matrix_list.size();
    size_t true_pos = 0;
    size_t false_pos = 0;
    size_t false_neg = 0;
    size_t true_neg = 0;
    size_t imports = 0;
    size_t unobsv_case = 0;

    Rcpp::NumericMatrix contacts;
    
    for (size_t m = 0; m < list_size; m++) {

      if (eps[m] < 0.0 || eta[m] < 0.0 || lambda[m] < 0.0) {
	return R_NegInf;
      }
	
      true_pos = false_pos = false_neg = true_neg = imports = unobsv_case = 0;
	  
      contacts = Rcpp::as<Rcpp::NumericMatrix>(ctd_matrix_list[m]);
      for (size_t j = 0; j < N; j++) {
	if (alpha[j] == NA_INTEGER) {
	  imports += 1;
	} else if (kappa[j] > 1) {
	  unobsv_case += 1;
	} else {
	  true_pos += contacts(j, alpha[j] - 1); // offset
	}
      }

      false_pos = C_nrow[m] - true_pos;
      false_neg = N - imports - unobsv_case - true_pos;
      true_neg = C_combn - true_pos - false_pos - false_neg;

      // untimed contact model
      if(lambda[m] == 0.0) {
	if(false_pos > 0) {
	  out += R_NegInf;
	} else {
	  out += log(eps[m]*eta[m]) * (double) true_pos +
	    log(1 - eps[m]*eta[m]) * (double) false_neg +
	    log(1 - eps[m]*lambda[m]) * (double) true_neg;
	}
      } else if(eta[m] == 0.0) {
	if(true_pos > 0) {
	  out += R_NegInf;
	} else {
	  out += log(eps[m]*lambda[m]) * false_pos +
	    log(1 - eps[m]*eta[m]) * false_neg +
	    log(1 - eps[m]*lambda[m]) * true_neg;
	}
      } else if(eta[m]*eps[m] == 1.0) {
	if(false_neg > 0) {
	  out += R_NegInf;
	} else {
	  out += log(eps[m]*eta[m]) * true_pos +
	    log(eps[m]*lambda[m]) * false_pos +
	    log(1 - eps[m]*lambda[m]) * true_neg;
	}
      } else if(lambda[m]*eps[m] == 1.0) {
	if(true_neg > 0) {
	  out += R_NegInf;
	} else {
	  out += log(eps[m]*eta[m]) * true_pos +
	    log(eps[m]*lambda[m]) * false_pos +
	    log(1 - eps[m]*eta[m]) * false_neg;
	}
      } else {
	out += log(eps[m]*eta[m]) * true_pos +
	  log(eps[m]*lambda[m]) * false_pos +
	  log(1 - eps[m]*eta[m]) * false_neg +
	  log(1 - eps[m]*lambda[m]) * true_neg;
	  
	// Rprintf("true pos %d | false pos %d | true neg %d | false neg %d\n",
	// 	  true_pos, false_pos, true_neg, false_neg);

      }
	
    } 
    
    return(out);
    
  } else { //use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(l[0]);
    int arity = l[1];
    if (arity == 3) return Rcpp::as<double>(f(data, param, i));
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

// This likelihood corresponds to the probability of observing a reported
// timeline (i.e. the place a given case is on a given day of the outbreak).

// The likelihood is based on the probability of a case in a given place i
// infecting a case in a given place j on the day of infection. 
//
// with:
// 'eps' is the probability of i = j
// 'tau' is probability of an unobserved case remaining in its place of infection

double cpp_ll_timeline(Rcpp::List data, Rcpp::List param, SEXP i,
		       Rcpp::RObject custom_function) {

  size_t N = static_cast<size_t>(data["N"]);
  if (N < 2) return 0.0;

  Rcpp::List l;
  if (custom_function != R_NilValue) {
    l = Rcpp::as<Rcpp::List>(custom_function);
  }
  if (custom_function == R_NilValue || (l.size() > 0 && l[0] == R_NilValue)) {

    Rcpp::List ctd_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_matrix"]);
    Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);

    Rcpp::List trans_mat_list = Rcpp::as<Rcpp::List>(param["trans_mat"]);
    
    if(ctd_timed_matrix_list.size() < 1) return 0.0;
    
    double out = 0;
    double p_wrong = data["p_wrong"];
    
    Rcpp::NumericVector eps = param["eps"];
    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector t_onw = param["t_onw"];
    
    Rcpp::IntegerVector N_place_vec = data["N_place"];

    Rcpp::NumericMatrix timeline;
    Rcpp::NumericVector trans_mat;


    int C_ind = static_cast<int>(data["C_ind"]);
    
    size_t list_size = ctd_timed_matrix_list.size();
    size_t imports = 0;
    size_t size_1;
    size_t size_2;
    int w1;
    int w2;
    int ind1;
    int ind2;
    size_t length_i;
    Rcpp::IntegerVector vec_i;
    size_t ind;
    int N_place;
    size_t j;

    if(i == R_NilValue) {
      // evaluate for all cases
      vec_i = Rcpp::seq(1, N);
      length_i = vec_i.size();
    } else {
      // only the cases listed in 'i' are retained
      // pass via tmp to so you don't run into declaration issues with vec_i
      Rcpp::IntegerVector tmp(i);
      vec_i = tmp;
      length_i = vec_i.size();
    }
    
    for (size_t m = 0; m < list_size; m++) {

      // index for eps
      ind = m + ctd_matrix_list.size();

      timeline = Rcpp::as<Rcpp::NumericMatrix>(ctd_timed_matrix_list[m]);
      trans_mat = Rcpp::as<Rcpp::NumericVector>(trans_mat_list[m]);

      N_place = N_place_vec[m] + 2;

      for (size_t k = 0; k < length_i; k++) {
	
	j = vec_i[k] - 1; // offset

	if (alpha[j] != NA_INTEGER) {

	  if(kappa[j] == 1) {
	    ind1 = ind2 = t_inf[j] + C_ind;
	  } else if(kappa[j] > 1) {
	    ind1 = t_inf[j] + C_ind;
	    ind2 = t_onw[j] + C_ind;
	  }
	  
	  if(ind1 >= 0 && ind1 < timeline.ncol() && ind2 >= 0 && ind2 < timeline.ncol()) {
	    w1 = static_cast<int>(timeline(j, ind1));
	    w2 = static_cast<int>(timeline(alpha[j] - 1, ind2));
	    if(w1 != -1 && w2 != -1) {
	      // if(j+1 == 3) {
	      // 	std::printf("YESCORRECT: t_inf = %i | ind1 = %i | i = %i | alpha_i = %i | w1 = %i | w2 = %i | ind = %i | ll = %f\n", t_inf[j], ind1, j+1, alpha[j], w1, w2, m, log(trans_mat(N_place*N_place*(kappa[j]-1) + N_place*(w1) + w2)));
	      // }
	      out += log(trans_mat(N_place*N_place*(kappa[j]-1) + N_place*(w1) + w2));
	    } else {
	      // if(j+1 == 3) {
	      // 	std::printf("INCORRECT: t_inf = %i | ind1 = %i | i = %i | alpha_i = %i | w1 = %i | w2 = %i | ind = %i | ll  %f\n", t_inf[j], ind1, j+1, alpha[j], w1, w2, m, log(p_wrong));
	      // }
	      out +=  log(p_wrong);
	    }
	  } else {
	    out +=  log(p_wrong);
	  }
	}
	
      }
    }
    
    return(out);
    
  } else { //use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(l[0]);
    int arity = l[1];
    if (arity == 3) return Rcpp::as<double>(f(data, param, i));
    return Rcpp::as<double>(f(data, param));
  }
  
}

double cpp_ll_timeline(Rcpp::List data, Rcpp::List param, size_t i,
		       Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timeline(data, param, si, custom_function);
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
      cpp_ll_timing_sampling(data, param, i) +
      cpp_ll_genetic(data, param, i) +
      cpp_ll_reporting(data, param, i);
  } else { // use of a customized likelihood functions
    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);
    return cpp_ll_timing_infections(data, param, i, list_functions["timing_infections"]) +
      cpp_ll_timing_sampling(data, param, i, list_functions["timing_sampling"]) +
      cpp_ll_genetic(data, param, i, list_functions["genetic"]) +
      cpp_ll_reporting(data, param, i, list_functions["reporting"]);
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
// - p(timline): see function cpp_ll_timeline

double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i,
		  Rcpp::RObject custom_functions) {

  if (custom_functions == R_NilValue) {

    return cpp_ll_timing_infections(data, param, i) +
      cpp_ll_timing_sampling(data, param, i) +
      cpp_ll_genetic(data, param, i) +
      cpp_ll_reporting(data, param, i) +
      cpp_ll_contact(data, param, i) +
      cpp_ll_timeline(data, param, i);

  }  else { // use of a customized likelihood functions
    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);

    return cpp_ll_timing_infections(data, param, i, list_functions["timing_infections"]) +
      cpp_ll_timing_sampling(data, param, i, list_functions["timing_sampling"]) +
      cpp_ll_genetic(data, param, i, list_functions["genetic"]) +
      cpp_ll_reporting(data, param, i, list_functions["reporting"]) +
      cpp_ll_contact(data, param, i, list_functions["contact"]) +
      cpp_ll_timeline(data, param, i, list_functions["timeline"]);

  }
}


double cpp_ll_all(Rcpp::List data, Rcpp::List param, size_t i,
              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_all(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}
