#include <Rcpp.h>
#include <Rmath.h>
#include "internals.h"
#include "likelihoods.h"


/*
  Movement of infection dates are +/- 1 from current states. These movements are currently
  vectorised, i.e. a bunch of dates are proposed all together; this may not be sustainable for
  larger datasets. The non-vectorised option will be slower and speed-up with C/C++ will be more
  substantial then.

  This version differs from the initial R implementation in several points:
  1. all cases are moved
  2. cases are moved one by one
  3. movement for each case is +/- 1 time unit

  Note that when computing the timing log-likelihood, the descendents of each case are also affected.

*/


// [[Rcpp::export("cpp.move.t.inf", rng = true)]]
void cpp_move_t_inf(Rcpp::List data, Rcpp::List param) {
  // Rcpp::List new_param(param);
  size_t N = static_cast<size_t>(data["N"]);
  size_t i = 0, j = 0, old_t_inf = 0, new_t_inf = 0;
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  Rcpp::IntegerVector t_inf = param["current.t.inf"]; // pointer to t_inf in param

  for (i = 0; i < N; i++) {
    // loglike with current value
    // old_loglike = cpp_ll_timing(data, param, i); // term for case 'i'
    // old_loglike = cpp_ll_timing(data, param, find_descendents(i)); // term descendents of 'i'
    old_loglike = cpp_ll_timing(data, param, R_NilValue);

    // proposal (+/- 1)
    old_t_inf = t_inf[i]; // save old value
    t_inf[i] += R::runif(0,1) > 0.5 ? 1 : -1; // new proposed value

   // loglike with current value
    // new_loglike = cpp_ll_timing(data, param, i); // term for case 'i'
    // new_loglike = cpp_ll_timing(data, param, find_descendents(i)); // term descendents of 'i'
    new_loglike = cpp_ll_timing(data, param, R_NilValue);

    // acceptance term
    p_accept = exp(new_loglike - old_loglike);

    // acceptance: the new value is already in t_inf, so we only act if the move is rejected, in
    // which case we restore the previous ('old') value
    if (p_accept < R::runif(0,1)) { // reject new values
      t_inf[i] = old_t_inf;
    }
  }

}
