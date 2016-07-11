#include <Rcpp.h>

/*
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
// [[Rcpp::export(rng = false)]]
double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i) {
  double mu = Rcpp::as<double>(param["current.mu"]);

  long int L = Rcpp::as<int>(data["L"]);

  Rcpp::NumericMatrix D = data["D"];
  Rcpp::IntegerVector alpha = param["current.alpha"];

  size_t length_nmut = 0, sum_nmut = 0;

  double * cD = REAL(D);
  size_t n = static_cast<size_t>(D.nrow());

  // all cases are retained
  if (i == R_NilValue) {
    for (size_t j = 0; j < n; j++) {
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += cD[j * n + alpha[j] - 1];
	length_nmut++;
      }
    }

  } else {
    // only the cases listed in 'i' are retained
    size_t ni = static_cast<size_t>(LENGTH(i));
    int * ci = INTEGER(i);
    for (size_t k = 0; k < ni; k++) {
      size_t j = ci[k] - 1;
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += cD[j * n + alpha[j] - 1];
	length_nmut++;
      }

    }
  }
  return log(mu / (1 - mu)) * sum_nmut + length_nmut * log(1 - mu) * L;
}
