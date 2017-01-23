#include <Rcpp.h>
#include <Rmath.h>
#include "internals.h"
#include "likelihoods.h"
#include "priors.h"



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

// Movement of the mutation rate 'mu' is done using a dumb normal proposal. This
// is satisfying for now - we only reject a few non-sensical values outside the
// range [0;1]. The SD of the proposal (implicitely contained in rand$mu.rnorm1,
// but really provided through 'config', seems fine as the range of real values
// will never change much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_mu(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.
  Rcpp::List new_param = clone(param); 
  Rcpp::NumericVector mu = param["mu"];
  Rcpp::NumericVector new_mu = new_param["mu"];

  double sd_mu = static_cast<double>(config["sd_mu"]);

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  // proposal (normal distribution with SD: config$sd_mu)

  new_mu[0] += R::rnorm(0.0, sd_mu); // new proposed value


  // automatic rejection of negative mu
  
  if (new_mu[0] < 0.0) {
    return param;
  }

   
  // compute likelihoods
  old_logpost = cpp_ll_genetic(data, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_genetic(data, new_param, R_NilValue, custom_ll);


  // compute priors

  old_logpost += cpp_prior_mu(param, config, custom_prior);
  new_logpost += cpp_prior_mu(new_param, config, custom_prior);
    

  // acceptance term

  p_accept = exp(new_logpost - old_logpost);


  // acceptance: the new value is already in mu, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value

  if (p_accept < unif_rand()) { // reject new values
    return param;
  }
  
  return new_param;
}






// ---------------------------

// Movement of the Reporting probability 'pi' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal (implicitely contained
// in rand$pi.rnorm1, but really provided through 'config', seems fine as the
// range of real values will never change much. Probably not much point in using
// auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_pi(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param); 
  Rcpp::NumericVector pi = param["pi"]; // these are just pointers
  Rcpp::NumericVector new_pi = new_param["pi"]; // these are just pointers

  double sd_pi = static_cast<double>(config["sd_pi"]);

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  
  // proposal (normal distribution with SD: config$sd_pi)
  
  new_pi[0] += R::rnorm(0.0, sd_pi); // new proposed value


  // automatic rejection of pi outside [0;1]

  if (new_pi[0] < 0.0 || new_pi[0] > 1.0) {
    return param;
  }

  
  // compute likelihoods
  old_logpost = cpp_ll_reporting(data, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_reporting(data, new_param, R_NilValue, custom_ll);


  // compute priors

  old_logpost += cpp_prior_pi(param, config, custom_prior);
  new_logpost += cpp_prior_pi(new_param, config, custom_prior);

  
  // acceptance term
  
  p_accept = exp(new_logpost - old_logpost);


  // acceptance: the new value is already in pi, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value
  
  if (p_accept < unif_rand()) { // reject new values
    return param;
  }
  
  return new_param;
}






// ---------------------------

// Movement of infection dates are +/- 1 from current states. These movements
// are currently vectorised, i.e. a bunch of dates are proposed all together;
// this may not be sustainable for larger datasets. The non-vectorised option
// will be slower and speed-up with C/C++ will be more substantial then.

// This version differs from the initial R implementation in several points:

// 1. all cases are moved
// 2. cases are moved one by one
// 3. movement for each case is +/- 1 time unit

// Notes

// - when computing the timing log-likelihood, the descendents of each
// case are also affected.

// - we generate a new vector 'new_t_inf', which will replace the
// previous pointer defining param["t_inf"].

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_t_inf(Rcpp::List data, Rcpp::List param,
			  Rcpp::RObject list_custom_ll = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.
  
  Rcpp::List new_param = clone(param); 
  Rcpp::IntegerVector t_inf = param["t_inf"];
  Rcpp::IntegerVector new_t_inf = new_param["t_inf"];
  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::IntegerVector local_cases;

  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
  double old_loc_loglike = 0.0, new_loc_loglike = 0.0, p_loc_accept = 0.0;

  
  for (size_t i = 0; i < N; i++) {
    // NOTE: local likelihood does not work here. Need to investigate why. 

    local_cases = cpp_find_descendents(param["alpha"], i+1);

    // loglike with current value
    // old_loglike = cpp_ll_timing(data, param, R_NilValue);
    old_loc_loglike = cpp_ll_timing(data, param, i+1, list_custom_ll); // term for case 'i' with offset

    // term descendents of 'i'
    if (local_cases.size() > 0) {
      old_loc_loglike += cpp_ll_timing(data, param, local_cases, list_custom_ll);
    }

    // proposal (+/- 1)
    new_t_inf[i] += unif_rand() > 0.5 ? 1 : -1; // new proposed value

    // loglike with new value
    // new_loglike = cpp_ll_timing(data, new_param, R_NilValue);
    new_loc_loglike = cpp_ll_timing(data, new_param, i+1, list_custom_ll); // term for case 'i' with offset

    // term descendents of 'i'
    if (local_cases.size() > 0) {
      new_loc_loglike += cpp_ll_timing(data, new_param, local_cases, list_custom_ll);
    }


    // acceptance term
    // p_accept = exp(new_loglike - old_loglike);
    p_loc_accept = exp(new_loc_loglike - old_loc_loglike);


    // acceptance: the new value is already in t_inf, so we only act if the move
    // is rejected, in which case we restore the previous ('old') value

    if (p_loc_accept < unif_rand()) { // reject new values
      new_t_inf[i] = t_inf[i];
    }
  }

  return new_param;

}






// ---------------------------

// Movement of ancestries ('alpha') is not vectorised, movements are made one
// case at a time. This procedure is simply about picking an infector at random
// amongst cases preceeding the case considered. The original version in
// 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't_inf', but
// current implementation is simpler and seems to mix at least as well. Proper
// movement of 'alpha' needs this procedure as well as a swapping procedure
// (swaps are not possible through move.alpha only); in all instances, 'alpha'
// is on the scale 1:N.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_alpha(Rcpp::List data, Rcpp::List param, 
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector new_alpha = new_param["alpha"];

  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;


  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved, if there is at least 1 choice
    if (alpha[i] != NA_INTEGER && sum(t_inf < t_inf[i]) > 0) {

      // loglike with current value
      // old_loglike = cpp_ll_all(data, param, R_NilValue);
      old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset

      // proposal (+/- 1)
      new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1); // new proposed value (on scale 1 ... N)

      // loglike with current value
      new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

      // acceptance term
      p_accept = exp(new_loglike - old_loglike);

      // which case we restore the previous ('old') value
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
      }
    }
  }

  return new_param;
}





// ---------------------------

// The basic movement of ancestries (picking an ancestor at random amongst in
// previous cases) makes swaps of ancestries (A->B) to (B->A) very
// difficult. This function addresses the issue. It is computer-intensive, but
// likely a determining factor for faster mixing. Unlike previous versions in
// the original 'outbreaker' package, all cases are 'moved' here. A swap is
// defined as:

// x -> a -> b  becomes a -> x -> b

// Obviously cases are moved one at a time. We need to use local likelihood
// changes for this move to scale well with outbreak size. The complicated bit
// is that the move impacts all descendents from 'a' as well as 'x'.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_swap_cases(Rcpp::List data, Rcpp::List param,
			       Rcpp::RObject list_custom_ll = R_NilValue) {

  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::List swapinfo; // contains alpha and t_inf
  Rcpp::IntegerVector local_cases;
  
  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;


  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved, if there is at least 1 choice
    if (alpha[i] != NA_INTEGER && sum(t_inf < t_inf[i]) > 0) {

      // The local likelihood is defined as the likelihood computed for the
      // cases affected by the swap; these include:

      // - 'i'
      // - the descendents of 'i'
      // - 'alpha[i]'
      // - the descendents of 'alpha[i]' (other than 'i')
      
      local_cases = cpp_find_local_cases(param["alpha"], i+1);
	

      // loglike with current parameters

      old_loglike = cpp_ll_all(data, param, local_cases, list_custom_ll); // offset


      // proposal: swap case 'i' and its ancestor

      swapinfo = cpp_swap_cases(param, i+1);
      new_param["alpha"] = swapinfo["alpha"];
      new_param["t_inf"] = swapinfo["t_inf"];
        
      
      // loglike with new parameters

      new_loglike = cpp_ll_all(data, new_param, local_cases, list_custom_ll);
      
      
      // acceptance term
      
      p_accept = exp(new_loglike - old_loglike);


      // acceptance: change param only if new values is accepted
      
      if (p_accept >= unif_rand()) { // accept new parameters
	param["alpha"] = new_param["alpha"];
	param["t_inf"] = new_param["t_inf"];
      }
    }
  }

  return param;
}







// ---------------------------


// Movement of the number of generations on transmission chains ('kappa') is
// done for one ancestry at a time. As for infection times ('t_inf') we use a
// dumb, symmetric +/- 1 proposal. But because values are typically in a short
// range (e.g. [1-3]) we probably propose more dumb values here. We may
// eventually want to bounce back or use and correct for assymetric proposals.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_kappa(Rcpp::List data, Rcpp::List param, Rcpp::List config, 
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector new_kappa = new_param["kappa"];

  size_t N = static_cast<size_t>(data["N"]);
  size_t K = static_cast<size_t>(config["max_kappa"]);
  size_t jump;

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved    
    if (alpha[i] != NA_INTEGER) {

      // propose new kappa
      jump = (unif_rand() > 0.5) ? 1 : -1;
      new_kappa[i] = kappa[i] + jump;


      // only look into this move if new kappa is positive and smaller than the
      // maximum value; if not, remember to reset the value of new_kappa to that
      // of kappa, otherwise we implicitely accept stupid moves automatically
      
      if (new_kappa[i] < 1 || new_kappa[i] > K) {
	new_kappa[i] = kappa[i];
      } else {
	
	// loglike with current parameters	
	old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll);


	// loglike with new parameters
	new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

				 
	// acceptance term
	p_accept = exp(new_loglike - old_loglike);


	// acceptance: change param only if new values is accepted      
	if (p_accept >= unif_rand()) { // accept new parameters
	  // Rprintf("\naccepting kappa:%d  (p: %f  old ll:  %f  new ll: %f",
	  // 		new_kappa[i], p_accept, old_loglike, new_loglike);
	  param["kappa"] = new_kappa;
	} else {
	  new_kappa[i] = kappa[i];
	}
      } 
    }

  }

  return param;
}
