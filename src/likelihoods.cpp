#include <Rcpp.h>

/*
  This likelihood corresponds to the probability of observing a number of mutations between cases
  and their ancestors. See src/likelihoods.cpp for details of the Rcpp implmentation.

  The likelihood is based on the number of mutations between a case and its ancestor;
  these are extracted from a pairwise genetic distance matrix (data$D)
  the log-likelihood is computed as: sum(mu^nmut + (1-mu)^(L-nmut))
  with:
  'mu' is the mutation probability
  'L' the number of sites in the alignment
  'nmut' the number of mutations between an ancestor and its descendent

  For computer efficiency, we re-factorise it as:
  log(mu / (1 - mu)) * sum(nmut) + length(nmut) * log(1 - mu) * L
  which limits to 2 operations rather than 2*n
  (tip from Rich Fitzjohn)

*/
// [[Rcpp::export("cpp.ll.genetic", rng = false)]]
double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i) {
  double mu = Rcpp::as<double>(param["current.mu"]);

  long int L = Rcpp::as<int>(data["L"]);

  Rcpp::NumericMatrix D = data["D"];
  Rcpp::IntegerVector alpha = param["current.alpha"]; // remember the '-1' offset!

  size_t length_nmut = 0, sum_nmut = 0;

  size_t N = static_cast<size_t>(data["N"]);

  // all cases are retained
  if (i == R_NilValue) {
    for (size_t j = 0; j < N; j++) {
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += D(j, alpha[j] - 1);
	length_nmut++;
      }
    }

  } else {
    // only the cases listed in 'i' are retained
    size_t length_i = static_cast<size_t>(LENGTH(i));
    Rcpp::IntegerVector vec_i(i);
    for (size_t k = 0; k < length_i; k++) {
      size_t j = vec_i[k] - 1;
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += D(j, alpha[j] - 1);
	length_nmut++;
      }

    }
  }
  return log(mu / (1 - mu)) * sum_nmut + length_nmut * log(1 - mu) * L;
}





/*

  This likelihood corresponds to the probability of observing infection dates of cases given the
  infection dates of their ancestors.
*/
// [[Rcpp::export("cpp.ll.timing.infections", rng = false)]]
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i) {
  size_t N = static_cast<size_t>(data["N"]);
  size_t j = 0;
  size_t delay = 0;

  Rcpp::IntegerVector alpha = param["current.alpha"]; // remember the '-1' offset!
  Rcpp::IntegerVector t_inf = param["current.t.inf"];
  Rcpp::IntegerVector kappa = param["current.kappa"];
  Rcpp::NumericMatrix w_dens = data["log.w.dens"];

  double out = 0.0;

  // all cases are retained
  if (i == R_NilValue) {
    for (j = 0; j < N; j++) {
      if (alpha[j] != NA_INTEGER) {
	delay = t_inf[j] - t_inf[alpha[j] - 1];
	if (delay < 1 || delay > w_dens.ncol()) {
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
      j = vec_i[k] - 1;
      if (alpha[j] != NA_INTEGER) {
	delay = t_inf[j] - t_inf[alpha[j] - 1];
	if (delay < 1 || delay > w_dens.ncol()) {
	  return  R_NegInf;
	}

	out += w_dens(kappa[j] - 1, delay - 1);
      }

    }
  }

  return out;
}






// This likelihood corresponds to the probability of reporting dates of cases given their
// infection dates.
// [[Rcpp::export("cpp.ll.timing.sampling", rng = false)]]
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i) {
  size_t N = static_cast<size_t>(data["N"]);
  size_t j = 0;
  size_t delay = 0;

  Rcpp::IntegerVector dates = data["dates"];
  Rcpp::IntegerVector t_inf = param["current.t.inf"];
  Rcpp::NumericVector f_dens = data["log.f.dens"];

  double out = 0.0;

  // all cases are retained
  if (i == R_NilValue) {
    for (j = 0; j < N; j++) {
      delay = dates[j] - t_inf[j];
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
      j = vec_i[k] - 1;
      delay = dates[j] - t_inf[j];
      if (delay < 1 || delay > f_dens.size()) {
	return  R_NegInf;
      }
      out += f_dens[delay - 1];
     }
  }

  return out;
}

