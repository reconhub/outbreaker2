#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <random> 
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
Rcpp::List cpp_move_mu(Rcpp::List param, Rcpp::List data, Rcpp::List config,
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

  // printf("old = %f | new = %f | ll_old = %f | ll_new = %f\n",
  // 	 mu[0], new_mu[0], old_logpost, new_logpost);

  // acceptance: the new value is already in mu, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value

  if (p_accept < unif_rand()) { // reject new values
    return param;
  }

  return new_param;
}






// ---------------------------

// movement of the Reporting probability 'pi' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal (implicitely contained
// in rand$pi.rnorm1, but really provided through 'config', seems fine as the
// range of real values will never change much. Probably not much point in using
// auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_pi(Rcpp::List param, Rcpp::List data, Rcpp::List config,
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

// movement of the Reporting probability 'pi' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal (implicitely contained
// in rand$pi.rnorm1, but really provided through 'config', seems fine as the
// range of real values will never change much. Probably not much point in using
// auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_tau(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			Rcpp::RObject custom_ll = R_NilValue,
			Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector tau = param["tau"]; // these are just pointers
  Rcpp::NumericVector eps = param["eps"]; // these are just pointers
  Rcpp::NumericVector new_tau = new_param["tau"]; // these are just pointers
  Rcpp::LogicalVector move_tau = config["move_tau"]; // these are just pointers
  Rcpp::NumericVector prop_u = data["prop_place_observed"]; // these are just pointers

  // transition matrics to be updated
  Rcpp::List new_trans_mat = Rcpp::as<Rcpp::List>(new_param["trans_mat"]);

  int max_kappa = static_cast<int>(config["max_kappa"]);
  double prop_tau_move = config["prop_tau_move"];

  // this is required to adjust the indexing when accessing epsilon
  int n_undated = eps.size() - tau.size();

  // raw transition probabilities
  Rcpp::List p_trans_list = Rcpp::as<Rcpp::List>(data["pp_trans_adj"]);
  Rcpp::List p_place_list = Rcpp::as<Rcpp::List>(data["pp_place"]);
  Rcpp::List p_place_adj_list = Rcpp::as<Rcpp::List>(data["pp_place_adj"]);

  Rcpp::NumericMatrix p_trans;
  Rcpp::NumericVector p_place;
  Rcpp::NumericVector p_place_adj;
  
  Rcpp::NumericVector sd_tau = config["sd_tau"];

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  for (size_t i = 0; i < tau.size(); i++) {

    if(move_tau[i] && unif_rand() < prop_tau_move) {

      // new proposed value
      new_tau[i] += R::rnorm(0.0, sd_tau[i]); 

      // automatic rejection of tau outside [0;1]
      if (new_tau[i] < 0.0 || new_tau[i] > 1.0) {
	return param;
      }

      p_trans = Rcpp::as<Rcpp::NumericMatrix>(p_trans_list[i]);
      p_place = Rcpp::as<Rcpp::NumericVector>(p_place_list[i]);
      p_place_adj = Rcpp::as<Rcpp::NumericVector>(p_place_adj_list[i]);

      // update trans mat
      new_trans_mat[i] = get_transition_mat(p_trans,
					    p_place,
					    p_place_adj,
					    eps[i + n_undated],
					    tau[i],
					    prop_u[i],
					    max_kappa);
      
      // compute likelihoods
      old_logpost = cpp_ll_timeline(data, param, R_NilValue, custom_ll);
      new_logpost = cpp_ll_timeline(data, new_param, R_NilValue, custom_ll);

      // compute priors
      old_logpost += cpp_prior_tau(param, config, custom_prior);
      new_logpost += cpp_prior_tau(new_param, config, custom_prior);


      // acceptance term
      p_accept = exp(new_logpost - old_logpost);

      // acceptance: the new value is already in tau, so we only act if the move is
      // rejected, in which case we restore the previous ('old') value

      if (p_accept < unif_rand()) { // reject new values
	new_tau[i] = tau[i];
      }
    
    }

  }
    
  return new_param;
  
}






// ---------------------------

// movement of the contact reporting coverage 'eps' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal is provided through
// 'config'; this seems fine as the range of real values will never change
// much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_eps(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			Rcpp::RObject list_custom_ll = R_NilValue,
			Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector eps = param["eps"]; // these are just pointers
  Rcpp::NumericVector tau = param["tau"]; // these are just pointers
  Rcpp::NumericVector new_eps = new_param["eps"]; // these are just pointers

  // transition matrics to be updated
  Rcpp::List new_trans_mat = Rcpp::as<Rcpp::List>(new_param["trans_mat"]);

  // old transition matrics
  Rcpp::List trans_mat = Rcpp::as<Rcpp::List>(param["trans_mat"]);
  
  // raw transition probabilities
  Rcpp::List p_trans_list = Rcpp::as<Rcpp::List>(data["pp_trans_adj"]);
  Rcpp::List p_place_list = Rcpp::as<Rcpp::List>(data["pp_place"]);
  Rcpp::List p_place_adj_list = Rcpp::as<Rcpp::List>(data["pp_place_adj"]);

  // get config
  Rcpp::LogicalVector move_eps = config["move_eps"];
  Rcpp::NumericVector sd_eps = config["sd_eps"];
  Rcpp::NumericVector prop_u = data["prop_place_observed"]; // these are just pointers
  int max_kappa = static_cast<int>(config["max_kappa"]);
  double prop_eps_move = config["prop_eps_move"];

  Rcpp::NumericMatrix p_trans;
  Rcpp::NumericVector p_place;
  Rcpp::NumericVector p_place_adj;
  
  // number of untimed contact types
  int n_ctd_untimed = eps.size() - p_trans_list.size();

  int ind;
  Rcpp::List list_func;
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  for (size_t i = 0; i < eps.size(); i++) {

    ind = i - n_ctd_untimed;
    
    if((move_eps[i] && ind < 0) ||
       (move_eps[i] && ind >= 0 && unif_rand() < prop_eps_move)) {
      
      // proposal (normal distribution with SD: config$sd_eps)

      new_eps[i] += R::rnorm(0.0, sd_eps[i]); // new proposed value

      // automatic rejection of eps outside [0;1]

      if (new_eps[i] < 0.0 || new_eps[i] > 1.0) {
	return param;
      }

      // if timed contacts available, update the transition matrices
      
      if(ind >= 0) {
	
	p_trans = Rcpp::as<Rcpp::NumericMatrix>(p_trans_list[ind]);
	p_place = Rcpp::as<Rcpp::NumericVector>(p_place_list[ind]);
	p_place_adj = Rcpp::as<Rcpp::NumericVector>(p_place_adj_list[ind]);
	
	new_trans_mat[ind] = get_transition_mat(p_trans,
						p_place,
						p_place_adj,
						new_eps[i],
						tau[ind],
						prop_u[ind],
						max_kappa);

	// use timeline likelihood
	if(list_custom_ll == R_NilValue) {
	  old_logpost = cpp_ll_timeline(data, param, R_NilValue);
	  new_logpost = cpp_ll_timeline(data, new_param, R_NilValue);
	} else {
	  // custom functions
	  list_func = Rcpp::as<Rcpp::List>(list_custom_ll);
	  old_logpost = cpp_ll_timeline(data, param, R_NilValue, list_func["timeline"]);
	  new_logpost = cpp_ll_timeline(data, new_param, R_NilValue, list_func["timeline"]);
	}
	
      } else {
	
	// use contacts likelihood
	if(list_custom_ll == R_NilValue) {
	  old_logpost = cpp_ll_contact(data, param, R_NilValue);
	  new_logpost = cpp_ll_contact(data, new_param, R_NilValue);
	} else {
	  list_func = Rcpp::as<Rcpp::List>(list_custom_ll);
	  old_logpost = cpp_ll_contact(data, param, R_NilValue, list_func["contact"]);
	  new_logpost = cpp_ll_contact(data, new_param, R_NilValue, list_func["contact"]);
	}
	
      }

      // compute priors
      old_logpost += cpp_prior_eps(param, config, custom_prior);
      new_logpost += cpp_prior_eps(new_param, config, custom_prior);

      // acceptance term
      p_accept = exp(new_logpost - old_logpost);

      // acceptance: the new value is already in eps, so we only act if the move is
      // rejected, in which case we restore the previous ('old') value

      if (p_accept < unif_rand()) { // reject new values
	// Rprintf("reject  %i | old = %f | new = %f | p_old = %f | p_new = %f\n",
	//  	i, eps[i], new_eps[i], old_logpost, new_logpost);
	new_eps[i] = eps[i];
	if(ind >= 0) {
	  new_trans_mat[ind] = trans_mat[ind];
	}
      } else {
	// replace values in param (needed because we use param to calculate
	// old_ll for eps[i+1];
	// Rprintf("accept  %i | old = %f | new = %f | p_old = %f | p_new = %f\n",
	//  	i, eps[i], new_eps[i], old_logpost, new_logpost);
	eps[i] = new_eps[i];
	if(ind >= 0) {
	  trans_mat[ind] = new_trans_mat[ind];
	}
      }
      
    }
  }
  
  return new_param;
  
}






// ---------------------------

// movement of the contact sensitivity 'eta' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal is provided through
// 'config'; this seems fine as the range of real values will never change
// much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_eta(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			Rcpp::RObject custom_ll = R_NilValue,
			Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector eta = param["eta"]; // these are just pointers
  Rcpp::NumericVector new_eta = new_param["eta"]; // these are just pointers
  Rcpp::LogicalVector move_eta = config["move_eta"]; // these are just pointers

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;
  
  //  double sd_eta = static_cast<double>(config["sd_eta"]);
  Rcpp::NumericVector sd_eta = config["sd_eta"];

  for (size_t i = 0; i < eta.size(); i++) {
    
    if(move_eta[i]) {
      
      // proposal (normal distribution with SD: config$sd_eta)

      new_eta[i] += R::rnorm(0.0, sd_eta[i]); // new proposed value

      // automatic rejection of eta outside [0;1]

      if (new_eta[i] < 0.0 || new_eta[i] > 1.0) {
	
	new_eta[i] = eta[i];
	
      } else {

	// compute likelihoods
	old_logpost = cpp_ll_contact(data, param, R_NilValue, custom_ll);
	new_logpost = cpp_ll_contact(data, new_param, R_NilValue, custom_ll);

	// compute priors
	old_logpost += cpp_prior_eta(param, config, custom_prior);
	new_logpost += cpp_prior_eta(new_param, config, custom_prior);

	// acceptance term
	p_accept = exp(new_logpost - old_logpost);

	//	std::printf("%i | eta_old = %f | eta_new = %f | ll_old = %f | ll_new = %f\n", i, eta[i], new_eta[i], old_logpost, new_logpost);

	
	// acceptance: the new value is already in eta, so we only act if the move is
	// rejected, in which case we restore the previous ('old') value

	if (p_accept < unif_rand()) { // reject new values
	  new_eta[i] = eta[i];
	} 
	
      }

    } else {

      new_eta[i] = eta[i];
	
    }

  }
  
  return new_param;
  
}






// ---------------------------

// movement of the non-infectious contact rate 'eps' is done using a dumb
// normal proposal. This is satisfying for now - we only reject a few
// non-sensical values outside the range [0;1]. The SD of the proposal is
// provided through 'config'; this seems fine as the range of real values will
// never change much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_lambda(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			   Rcpp::RObject custom_ll = R_NilValue,
			   Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector lambda = param["lambda"]; // these are just pointers
  Rcpp::NumericVector new_lambda = new_param["lambda"]; // these are just pointers
  Rcpp::LogicalVector move_lambda = config["move_lambda"]; // these are just pointers

  //  double sd_lambda = static_cast<double>(config["sd_lambda"]);
  Rcpp::NumericVector sd_lambda = config["sd_lambda"];
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;


  for (size_t i = 0; i < lambda.size(); i++) {

    if(move_lambda[i]) {
      
      // proposal (normal distribution with SD: config$sd_lambda)

      new_lambda[i] += R::rnorm(0.0, sd_lambda[i]); // new proposed value


      // automatic rejection of lambda outside [0;1]

      if (new_lambda[i] < 0.0 || new_lambda[i] > 1.0) {
	return param;
      }


      // compute likelihoods
      old_logpost = cpp_ll_contact(data, param, R_NilValue, custom_ll);
      new_logpost = cpp_ll_contact(data, new_param, R_NilValue, custom_ll);


      // compute priors

      old_logpost += cpp_prior_lambda(param, config, custom_prior);
      new_logpost += cpp_prior_lambda(new_param, config, custom_prior);


      // acceptance term

      p_accept = exp(new_logpost - old_logpost);


      // acceptance: the new value is already in lambda, so we only act if the move is
      // rejected, in which case we restore the previous ('old') value

      if (p_accept < unif_rand()) { // reject new values
	new_lambda[i] = lambda[i];
      }

    } else {

      new_lambda[i] = lambda[i];
      
    }
  
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
Rcpp::List cpp_move_t_inf(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector t_inf = param["t_inf"];
  Rcpp::IntegerVector kappa = param["kappa"];
  Rcpp::IntegerVector new_t_inf = new_param["t_inf"];
  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::IntegerVector local_cases;
  size_t N = static_cast<size_t>(data["N"]);

  // standard deviation of normal proposal for infection times
  double sd_t_inf = static_cast<double>(config["sd_t_inf"]);

  // define the vectors to sample from  
  Rcpp::IntegerVector sample_from = Rcpp::seq(1, 50);

  // change in infection time
  int t_inf_change;
  
  // define probabilities
  Rcpp::NumericVector sample_from_prob = Rcpp::dnorm(sample_from, 0.0, sd_t_inf);
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
  double old_loc_loglike = 0.0, new_loc_loglike = 0.0, p_loc_accept = 0.0;

  Rcpp::List list_functions = Rcpp::as<Rcpp::List>(list_custom_ll);

  for (size_t i = 0; i < N; i++) {

    local_cases = cpp_find_descendents(param["alpha"], i+1);

    // loglike with current value
    // term for case 'i' with offset
    old_loc_loglike = cpp_ll_timing(data, param, i+1, list_custom_ll) +
      cpp_ll_timeline(data, param, i+1, list_functions["timeline"]);

    // term descendents of 'i'
    if (local_cases.size() > 0) {
      old_loc_loglike += cpp_ll_timing(data, param, local_cases, list_custom_ll);
    }
    
    // draw change in t_inf
    t_inf_change = Rcpp::as<int>(Rcpp::sample(sample_from, 1, true, sample_from_prob));

    // new proposed value (either add or subtract)
    new_t_inf[i] += unif_rand() > 0.5 ? t_inf_change : -1*t_inf_change;
    
    // loglike with new value
    // term for case 'i' with offset
    new_loc_loglike = cpp_ll_timing(data, new_param, i+1, list_custom_ll) +
      cpp_ll_timeline(data, new_param, i+1, list_functions["timeline"]);

    // term descendents of 'i'
    // t_inf moves don't affect local cases in ll_timeline
    if (local_cases.size() > 0) {
      new_loc_loglike += cpp_ll_timing(data, new_param, local_cases, list_custom_ll);
    }

    // acceptance term
    p_loc_accept = exp(new_loc_loglike - old_loc_loglike);
    
    // acceptance: the new value is already in t_inf, so we only act if the move
    // is rejected, in which case we restore the previous ('old') value
    if (p_loc_accept < unif_rand()) { // reject new values
      new_t_inf[i] = t_inf[i];
    } else {
      t_inf[i] = new_t_inf[i];
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
Rcpp::List cpp_move_alpha(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector new_alpha = new_param["alpha"];
  Rcpp::LogicalVector move_alpha = config["move_alpha"];

  // Genetic sequence model
  std::string genetic_model = Rcpp::as<std::string>(data["genetic_model"]);

  Rcpp::NumericMatrix ancestors = param["ancestors"];
  Rcpp::NumericMatrix new_ancestors = new_param["ancestors"];

  Rcpp::NumericMatrix mrca = param["mrca"];
  Rcpp::NumericMatrix new_mrca = new_param["mrca"];

  Rcpp::NumericMatrix combn = data["dna_combn"];
  Rcpp::IntegerVector has_dna = data["has_dna_ind"];
  
  // Deep copy to prevent pointer exchange
  
  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
  double old_loc_loglike = 0.0, new_loc_loglike = 0.0;

  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved when move_alpha[i] is TRUE, if at least 1 choice
    if (alpha[i] != NA_INTEGER && move_alpha[i] && sum(t_inf < t_inf[i]) > 0) {

      // loglike with current value
      old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset

      // proposal (+/- 1)
      // new proposed value (on scale 1 ... N)
      new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1);

      if(genetic_model == "mrca") {
      
	// update ancestry matrix
	new_ancestors = cpp_find_ancestors(new_alpha, ancestors, has_dna);
      
	// // update mrca
	new_mrca = update_mrca(combn, new_ancestors);

	new_param["ancestors"] = new_ancestors;
	new_param["mrca"] = new_mrca;

      }
      
      // loglike with current value
      new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

      // acceptance term
      p_accept = exp(new_loglike - old_loglike);
      
      // which case we restore the previous ('old') value
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
	if(genetic_model == "mrca") {
	  new_ancestors = clone(ancestors);
	  new_mrca = clone(mrca);
	}
      } else {
	alpha[i] = new_alpha[i];
	if(genetic_model == "mrca") {
	  ancestors = clone(new_ancestors);
	  mrca = clone(new_mrca);
	}
      }      
    }
  }

  if(genetic_model == "mrca") {
    new_param["ancestors"] = ancestors;
    new_param["mrca"] = mrca;
  }
  
  return new_param;
  
}





// ---------------------------

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_model(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  
  Rcpp::List new_param = clone(param);

  // Genetic sequence model
  std::string genetic_model = Rcpp::as<std::string>(data["genetic_model"]);

  // Only propose a model swap prop_model_move % of the time
  double prop_model_move = config["prop_model_move"];
  Rcpp::LogicalVector move_alpha = config["move_alpha"];
  
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector t_onw = param["t_onw"]; // pointer to param$t_onw
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa

  size_t N = static_cast<size_t>(data["N"]);
  size_t K = static_cast<size_t>(config["max_kappa"]);
  
  Rcpp::IntegerVector new_t_onw = new_param["t_onw"];
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  Rcpp::IntegerVector new_alpha = new_param["alpha"];

  Rcpp::NumericMatrix ancestors = param["ancestors"];
  Rcpp::NumericMatrix new_ancestors = new_param["ancestors"];

  Rcpp::NumericMatrix mrca = param["mrca"];
  Rcpp::NumericMatrix new_mrca = new_param["mrca"];

  Rcpp::NumericMatrix combn = data["dna_combn"];
  Rcpp::IntegerVector has_dna = data["has_dna_ind"];

  // N_times is the number of days for which we have timeline data
  int N_times = data["N_times"];

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  int C_ind = static_cast<int>(data["C_ind"]);
  
  // N_place is the number of unique places for each type of timeline data
  Rcpp::IntegerVector N_place = data["N_place"];
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  bool m1_to_m2;
  Rcpp::IntegerVector t_seq;
  Rcpp::IntegerVector kappa_seq;
  int ind1;
  double tdiff;
  
  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved when move_alpha[i] is TRUE, if at least 1 choice
    if (alpha[i] != NA_INTEGER && move_alpha[i] &&
	sum(t_inf < t_inf[i]) > 0 &&
	prop_model_move < unif_rand()) {

      // loglike with current value
      // old_loglike = cpp_ll_all(data, param, R_NilValue);
      old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset

      // Pick a new sampled ancestor (new proposed value (on scale 1 ... N)
      new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1);

      if(genetic_model == "mrca") {

	// update ancestry matrix
	new_ancestors = cpp_find_ancestors(new_alpha, ancestors, has_dna);

	// update mrca
	new_mrca = update_mrca(combn, new_ancestors);

	new_param["ancestors"] = new_ancestors;
	new_param["mrca"] = new_mrca;

      }

      // Make sure m1 <-> m2 are proposed 50% of the time (otherwise you will
      // propose m1 -> m2 more if you spend more time in m1, and vice versa)
      m1_to_m2 = (unif_rand() > 0.5) ? true : false;

      // jump from Model1 to Model2
      if(kappa[i] == 1 && m1_to_m2) {

	// In Model2 (kappa > 1), infection times must be at least two days apart
	if(t_inf[i] - t_inf[new_alpha[i] - 1] < 2) {
	  p_accept = 0.0;
	} else {

	  // Propose t_onw between the infection times of the two sampled cases
	  t_seq = Rcpp::seq(t_inf[new_alpha[i] - 1] + 1, t_inf[i] - 1);
	  new_t_onw[i] = Rcpp::as<int>(Rcpp::sample(t_seq, 1, true));

	  ind1 = new_t_onw[i] + C_ind;
	  
	  // Propose kappa between the 2 and max_kappa
	  kappa_seq = Rcpp::seq(2, K);
	  new_kappa[i] = Rcpp::as<int>(Rcpp::sample(kappa_seq, 1, true));

	  // Propose a new place using prior probabilities. This ensures that
	  // our proposal probabilities and prior probabilites cancel out in the
	  // between model move. We also propose place = N_place + 1, which
	  // represents an unsampled place.
	  
	  // loglike with current value
	  new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);
	  
	  // The ratio of proposal distributions is length((t_inf_i +
	  // 1):(t_inf(alpha_i)-1)) The denominator is the uniform
	  // prior i.e. 1/N_times
	  if(new_loglike == R_NegInf) {
	    p_accept = 0.0;
	  } else {
	    // p_accept = exp(new_loglike - old_loglike)*t_seq.size()/N_times;
	    p_accept = exp(new_loglike - old_loglike);
	  }
	}
	// Jump from Model2 to Model1
      } else if(kappa[i] > 1 && !m1_to_m2) {

	// No longer inferring t_onw - set to -1000 (NA_INTEGER causes segfault)
	new_t_onw[i] = -1000;
	new_kappa[i] = 1;
	
	// loglike with current value
	new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);
	
	// if the new likelihood state is impossible, reject it - otherwise it
	// might get accepted because the ratio of proposal distribution says so 
	if(new_loglike == R_NegInf) {
	  p_accept = 0.0;
	} else {
	// The proposal probability of going from state M2 -> M1 is 1 (there is
	// only possible Model 1, in that we remove t_onw and kappa_beta_i); the
	// proposal probability of going from M1 -> M2 is 1/(range of possible
	// t_onw * range of possible kappa_beta_i) - we therefore need adjust
	// for this - however the prior on these infection times is 1/(Tend -
	// T0), which goes on the numerator
	  tdiff = t_inf[i] - t_inf[new_alpha[i] - 1] - 1;
	  if(tdiff == 0.0) {
	    p_accept = 0;
	  } else {
	    // p_accept = exp(new_loglike - old_loglike)*N_times/tdiff;
	    p_accept = exp(new_loglike - old_loglike);
	  }
	}
      } else {
	p_accept = 0.0;
      }

      // which case we restore the previous ('old') value
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
	new_t_onw[i] = t_onw[i];
	new_kappa[i] = kappa[i];
	if(genetic_model == "mrca") {
	  new_ancestors = clone(ancestors);
	  new_mrca = clone(mrca);
	}
      } else {
	alpha[i] = new_alpha[i];
	t_onw[i] = new_t_onw[i];
	kappa[i] = new_kappa[i];
	if(genetic_model == "mrca") {
	  mrca = clone(new_mrca);
	  ancestors = clone(new_ancestors);
	}
      }
    }
  }

  if(genetic_model == "mrca") {
    new_param["ancestors"] = ancestors;
    new_param["mrca"] = mrca;
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
Rcpp::List cpp_move_joint(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);

  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector t_onw = param["t_onw"];
  Rcpp::IntegerVector kappa = param["kappa"];

  // Genetic sequence model
  std::string genetic_model = Rcpp::as<std::string>(data["genetic_model"]);

  Rcpp::NumericMatrix ancestors = param["ancestors"];
  Rcpp::NumericMatrix new_ancestors = new_param["ancestors"];

  Rcpp::NumericMatrix mrca = param["mrca"];
  Rcpp::NumericMatrix new_mrca = new_param["mrca"];

  Rcpp::NumericMatrix combn = data["dna_combn"];
  Rcpp::IntegerVector has_dna = data["has_dna_ind"];
  
  size_t N = static_cast<size_t>(data["N"]);
  size_t K = static_cast<size_t>(config["max_kappa"]);
  Rcpp::LogicalVector move_alpha = config["move_alpha"];

  Rcpp::IntegerVector N_place = data["N_place"];
  double prop_alpha_move = config["prop_alpha_move"];

  int jump;

  double sd_t_onw = static_cast<double>(config["sd_t_onw"]);

  Rcpp::IntegerVector new_t_onw = new_param["t_onw"];
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  Rcpp::IntegerVector new_alpha = new_param["alpha"];
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
    
  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved when move_alpha[i] is TRUE, if at least 1 choice
    if (alpha[i] != NA_INTEGER && move_alpha[i] && sum(t_inf < t_inf[i]) > 0) {

      // loglike with current value

      // Within Model 1 move; just an alpha move (t_onw and kappa are not inferred)
      if(kappa[i] == 1) {

	old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset
	
	// proposal (+/- 1)
	// new proposed value (on scale 1 ... N)
	new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1); 

	if(genetic_model == "mrca") {

	  // update ancestry matrix
	  new_ancestors = cpp_find_ancestors(new_alpha, ancestors, has_dna);

	  // update mrca
	  new_mrca = update_mrca(combn, new_ancestors);

	  new_param["ancestors"] = new_ancestors;
	  new_param["mrca"] = new_mrca;

	}
	
	// loglike with current value
	new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

	// acceptance term
	p_accept = exp(new_loglike - old_loglike);

	// Within Model 2 move; move alpha, kappa and t_onw
      } else if(kappa[i] > 1) {

	old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset
	
	jump = (unif_rand() > 0.5) ? 1 : -1;
	new_kappa[i] = new_kappa[i] + jump;
	
	// A move back to kappa = 1 is impossible; that would represent a move to Model 1
	if (new_kappa[i] < 2 || new_kappa[i] > K) {
	  p_accept = 0.0;
	} else {
	  
	  // Only propose a new ancestry prop_alpha_move % of the time - this
	  // allows t_onw to mix independently of alpha at times
	  if(unif_rand() < prop_alpha_move) {
	    new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1);

	    if(genetic_model == "mrca") {
	      
	      // update ancestry matrix
	      new_ancestors = cpp_find_ancestors(new_alpha, ancestors, has_dna);

	      // update mrca
	      new_mrca = update_mrca(combn, new_ancestors);

	      new_param["ancestors"] = new_ancestors;
	      new_param["mrca"] = new_mrca;

	    }
	    
	  }
	  
	  // Normal MH move from previous state
	  new_t_onw[i] += std::round(R::rnorm(0.0, sd_t_onw));
	  
	  //	  loglike with current value
	  new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

	  // acceptance term
	  p_accept = exp(new_loglike - old_loglike);
	  
	}
      }

      // Reject new values if not accepted
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
	new_t_onw[i] = t_onw[i];
	new_kappa[i] = kappa[i];
	if(genetic_model == "mrca") {
	  new_ancestors = clone(ancestors);
	  new_mrca = clone(mrca);
	}
      } else {  // update param for next iteration if accepted 
	alpha[i] = new_alpha[i];
	t_onw[i] = new_t_onw[i];
	kappa[i] = new_kappa[i];
	if(genetic_model == "mrca") {
	  ancestors = clone(new_ancestors);
	  mrca = clone(new_mrca);
	}
      }
    }
  }
  
  if(genetic_model == "mrca") {
    new_param["ancestors"] = ancestors;
    new_param["mrca"] = mrca;
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
Rcpp::List cpp_move_swap_cases(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			       Rcpp::RObject list_custom_ll = R_NilValue) {

  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector t_onw = param["t_onw"]; // pointer to param$t_inf

  // Genetic sequence model
  std::string genetic_model = Rcpp::as<std::string>(data["genetic_model"]);
  
  Rcpp::NumericMatrix ancestors = param["ancestors"];
  Rcpp::NumericMatrix new_ancestors = new_param["ancestors"];

  Rcpp::NumericMatrix mrca = param["mrca"];
  Rcpp::NumericMatrix new_mrca = new_param["mrca"];

  Rcpp::NumericMatrix combn = data["dna_combn"];
  Rcpp::IntegerVector has_dna = data["has_dna_ind"];
  
  Rcpp::List swapinfo; // contains alpha and t_inf
  Rcpp::IntegerVector local_cases;

  bool swap_place = Rcpp::as<bool>(data["swap_place"]);
  size_t N = static_cast<size_t>(data["N"]);

  double prop_alpha_move = config["prop_alpha_move"];
  Rcpp::LogicalVector move_alpha = config["move_alpha"];

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  size_t n_mat = ctd_timed_matrix_list.size();

  // Shuffle indices to make equal cases equally likely
  Rcpp::IntegerVector idx = Rcpp::sample(N, N, false) - 1;

  for (size_t j = 0; j < N; j++) {

    size_t i = (size_t)idx[j];
    
    // only non-NA ancestries are moved when move_alpha[i] is TRUE, if at least 1 choice
    if (alpha[i] != NA_INTEGER && move_alpha[i] &&
	sum(t_inf < t_inf[i]) > 0 &&
	unif_rand() < prop_alpha_move) {

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

      swapinfo = cpp_swap_cases(param, i+1, swap_place);
      new_param["alpha"] = swapinfo["alpha"];
      new_param["t_inf"] = swapinfo["t_inf"];
      new_param["t_onw"] = swapinfo["t_onw"];
      new_param["kappa"] = swapinfo["kappa"];
      
      if(genetic_model == "mrca") {
	
	// update ancestry matrix
	new_ancestors = cpp_find_ancestors(new_param["alpha"], ancestors, has_dna);

	// update mrca
	new_mrca = update_mrca(combn, new_ancestors);

	new_param["ancestors"] = new_ancestors;
	new_param["mrca"] = new_mrca;

      }
      
      // loglike with new parameters
      new_loglike = cpp_ll_all(data, new_param, local_cases, list_custom_ll);


      // acceptance term

      p_accept = exp(new_loglike - old_loglike);
	          
      // acceptance: change param only if new values is accepted

      if (p_accept >= unif_rand()) { // accept new parameters
	param["alpha"] = new_param["alpha"];
	param["t_inf"] = new_param["t_inf"];
	param["t_onw"] = new_param["t_onw"];
	param["kappa"] = new_param["kappa"];
	if(genetic_model == "mrca") {
	  param["ancestors"] = clone(new_ancestors);
	  param["mrca"] = clone(new_mrca);
	  ancestors = clone(new_ancestors);
	  mrca = clone(new_mrca);
	}
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
Rcpp::List cpp_move_kappa(Rcpp::List param, Rcpp::List data, Rcpp::List config,
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
	  // new_kappa[i], p_accept, old_loglike, new_loglike);
	  kappa[i] = new_kappa[i];
	} else {
	  new_kappa[i] = kappa[i];
	}
      }
    }

  }

  return param;
}
