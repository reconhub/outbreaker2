#include <Rcpp.h>
#include <Rmath.h>


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
// 'nmut' the number of mutations between an ancestor and its descendent

// For computer efficiency, we re-factorise it as:
// log(mu / (1 - mu)) * sum(nmut) + length(nmut) * log(1 - mu) * L
// which limits to 2 operations rather than 2*n
// (tip from Rich Fitzjohn)

// [[Rcpp::export(rng = false)]]
double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i) {
  Rcpp::NumericMatrix D = data["D"];
  if (D.ncol() < 1) return 0.0;

  size_t N = static_cast<size_t>(data["N"]);
  if (N < 2) return 0.0;
  
  Rcpp::NumericMatrix w_dens = data["log.w.dens"];
  size_t K = w_dens.nrow();

  double mu = Rcpp::as<double>(param["mu"]);
  long int L = Rcpp::as<int>(data["L"]);
  Rcpp::IntegerVector alpha = param["alpha"]; // values are on 1:N
  Rcpp::IntegerVector kappa = param["kappa"];

  size_t length_nmut = 0, sum_nmut = 0;

  // p(mu < 0) = 0
  if (mu < 0.0) {
    return R_NegInf;
  }
  
  // all cases are retained
  if (i == R_NilValue) {
    for (size_t j = 0; j < N; j++) {
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += D(j, alpha[j] - 1); // offset
	length_nmut++;
      }
    }

  } else {
    // only the cases listed in 'i' are retained
    size_t length_i = static_cast<size_t>(LENGTH(i));
    Rcpp::IntegerVector vec_i(i);
    for (size_t k = 0; k < length_i; k++) {
      size_t j = vec_i[k] - 1; // offset
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += D(j, alpha[j] - 1); // offset
	length_nmut++;
      }

    }
  }
  return log(mu / (1 - mu)) * sum_nmut + length_nmut * log(1 - mu) * L;
}


double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, size_t i) {
  return cpp_ll_genetic(data, param, Rcpp::wrap(i));
}







// ---------------------------

// This likelihood corresponds to the probability of observing infection dates
// of cases given the infection dates of their ancestors.

// [[Rcpp::export(rng = false)]]
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::IntegerVector t_inf = param["t.inf"];
  Rcpp::IntegerVector kappa = param["kappa"];
  Rcpp::NumericMatrix w_dens = data["log.w.dens"];
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
}


double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, size_t i) {
  return cpp_ll_timing_infections(data, param, Rcpp::wrap(i));
}





// ---------------------------

// This likelihood corresponds to the probability of reporting dates of cases
// given their infection dates.

// [[Rcpp::export(rng = false)]]
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  Rcpp::IntegerVector dates = data["dates"];
  Rcpp::IntegerVector t_inf = param["t.inf"];
  Rcpp::NumericVector f_dens = data["log.f.dens"];

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
}


double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, size_t i) {
  return cpp_ll_timing_sampling(data, param, Rcpp::wrap(i));
}






// ---------------------------

// This likelihood corresponds to the probability of a given number of
// unreported cases on an ancestry.

// The likelihood is given by a geometric distribution with probability 'pi'
// to report a case

// - 'kappa' is the number of generation between two successive cases
// - 'kappa-1' is the number of unreported cases

// [[Rcpp::export(rng = false)]]
double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i) {
  Rcpp::NumericMatrix w_dens = data["log.w.dens"];
  size_t K = w_dens.nrow();
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;

  double pi = static_cast<double>(param["pi"]);
  Rcpp::IntegerVector kappa = param["kappa"];

  double out = 0.0;

  // p(pi < 0) = p(pi > 1) = 0
  if (pi < 0.0 || pi > 1.0) {
    return R_NegInf;
  }

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
}


double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, size_t i) {
  return cpp_ll_reporting(data, param, Rcpp::wrap(i));
}







// ---------------------------

// This likelihood corresponds to the sums of the separate timing likelihoods,
// which include:

// - p(infection dates): see function cpp_ll_timing_infections
// - p(collection dates): see function cpp_ll_timing_sampling

// [[Rcpp::export(rng = false)]]
double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i) {

  return cpp_ll_timing_infections(data, param, i) + 
    cpp_ll_timing_sampling(data, param, i);
}

double cpp_ll_timing(Rcpp::List data, Rcpp::List param, size_t i) {
  return cpp_ll_timing(data, param, Rcpp::wrap(i));
}





// ---------------------------

// This likelihood corresponds to the sums of the separate likelihoods, which
// include:

// - p(infection dates): see function cpp_ll_timing_infections
// - p(collection dates): see function cpp_ll_timing_sampling
// - p(genetic diversity): see function cpp_ll_genetic
// - p(missing cases): see function cpp_ll_reporting

// [[Rcpp::export(rng = false)]]
double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i) {

  return cpp_ll_timing_infections(data, param, i) +
    cpp_ll_timing_sampling(data, param, i) +
    cpp_ll_genetic(data, param, i) +
    cpp_ll_reporting(data, param, i);
}


double cpp_ll_all(Rcpp::List data, Rcpp::List param, size_t i) {
  return cpp_ll_all(data, param, Rcpp::wrap(i));
}
