#include <Rcpp.h>
#include <Rmath.h>
#include "internals.h"
#include "likelihoods.h"




/*

  Movement of the mutation rate 'mu' is done using a dumb normal proposal. This is satisfying for
  now - we only reject a few non-sensical values outside the range [0;1]. The SD of the proposal
  (implicitely contained in rand$mu.rnorm1, but really provided through 'config', seems fine as the
  range of real values will never change much. Probably not much point in using auto-tuning here.

*/

// [[Rcpp::export("cpp.move.mu", rng = true)]]
void cpp_move_mu(Rcpp::List data, Rcpp::List param, Rcpp::List config) {
  Rcpp::NumericVector mu = param["current.mu"]; // pointer to param$current.mu
  Rcpp::NumericVector new_mu = clone(mu); // new vector

  double sd_mu = static_cast<double>(config["sd.mu"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  // loglike with current value
  old_loglike = cpp_ll_genetic(data, param, R_NilValue);
  // Rcpp::Rcout << "old_mu:" << mu[0] << std::endl;
  // Rcpp::Rcout << "old_loglike:" << old_loglike << std::endl;

  // proposal (normal distribution with SD: config$sd.mu)
  new_mu[0] += R::rnorm(0, sd_mu); // new proposed value

  // loglike with current value
  param["current.mu"] = new_mu;
  new_loglike = cpp_ll_genetic(data, param, R_NilValue);
  // Rcpp::Rcout << "new_mu:" << new_mu[0] << std::endl;
  // Rcpp::Rcout << "new_loglike:" << new_loglike << std::endl;

  // acceptance term
  p_accept = exp(new_loglike - old_loglike);

  // Rcpp::Rcout << "p_accept:" << p_accept << std::endl;

  // acceptance: the new value is already in mu, so we only act if the move is rejected, in
  // which case we restore the previous ('old') value
  if (p_accept < R::runif(0,1)) { // reject new values
    param["current.mu"] = mu;
  }
}




/*

  Movement of infection dates are +/- 1 from current states. These movements are currently
  vectorised, i.e. a bunch of dates are proposed all together; this may not be sustainable for
  larger datasets. The non-vectorised option will be slower and speed-up with C/C++ will be more
  substantial then.

  This version differs from the initial R implementation in several points:
  1. all cases are moved
  2. cases are moved one by one
  3. movement for each case is +/- 1 time unit

  Notes
  - when computing the timing log-likelihood, the descendents of each case are also affected.
  - we generate a new vector 'new_t_inf', which will replace the previous pointer defining 
  param["current.t.inf"].

*/

// [[Rcpp::export("cpp.move.t.inf", rng = true)]]
void cpp_move_t_inf(Rcpp::List data, Rcpp::List param) {
  Rcpp::IntegerVector t_inf = param["current.t.inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector new_t_inf = clone(t_inf); // new vector

  size_t N = static_cast<size_t>(data["N"]);

  size_t i = 0, j = 0, old_t_inf_i = 0;

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;


  for (i = 0; i < N; i++) {
    // loglike with current value
    // old_loglike = cpp_ll_timing(data, param, i); // term for case 'i'
    // old_loglike = cpp_ll_timing(data, param, find_descendents(i)); // term descendents of 'i'
    old_loglike = cpp_ll_timing(data, param, R_NilValue);

    // proposal (+/- 1)
    new_t_inf[i] += R::runif(0,1) > 0.5 ? 1 : -1; // new proposed value

    // loglike with current value
    // new_loglike = cpp_ll_timing(data, param, i); // term for case 'i'
    // new_loglike = cpp_ll_timing(data, param, find_descendents(i)); // term descendents of 'i'
    param["current.t.inf"] = new_t_inf;
    new_loglike = cpp_ll_timing(data, param, R_NilValue);

    // acceptance term
    p_accept = exp(new_loglike - old_loglike);

    // acceptance: the new value is already in t_inf, so we only act if the move is rejected, in
    // which case we restore the previous ('old') value
    if (p_accept < R::runif(0,1)) { // reject new values
      new_t_inf[i] = t_inf[i];
      param["current.t.inf"] = t_inf;
    }
  }

  // as we haven't touched the content of t_inf, and all new values are in new_t_inf, we need to
  // make sure this new vector replaces the previous one
  param["current.t.inf"] = new_t_inf;
}


