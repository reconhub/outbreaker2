#include <Rcpp.h>

// [[Rcpp::export]]
double ll_genetic(Rcpp::List data, Rcpp::List param, Rcpp::IntegerVector i) {
  double mu = Rcpp::as<double>(data["mu"]);

  long int L = Rcpp::as<int>(data["L"]);

  Rcpp::NumericMatrix D = data["D"];
  Rcpp::IntegerVector alpha = param["current.alpha"];

  size_t length_nmut = 0, sum_nmut = 0;

  int * cD = INTEGER(D);
  size_t n = static_cast<size_t>(D.nrow());

  if (i.size() == 0) {
    for (size_t j = 0; j < n; j++) {
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += cD[j * n + alpha[j] - 1];
	length_nmut++;
      }
    }

  } else {
    size_t n = static_cast<size_t>(i.size());
    int * ci = INTEGER(i);
    for (size_t k = 0; k < n; k++) {
      size_t j = ci[k];
      if (alpha[j] != NA_INTEGER) {
	sum_nmut += cD[j * n + alpha[j] - 1];
	length_nmut++;
      }

    }
  }
  return log(mu / (1 - mu)) * sum_nmut + length_nmut * log(1 - mu) * L;
}
