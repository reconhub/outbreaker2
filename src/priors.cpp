#include <Rcpp.h>
#include <Rmath.h>
#include "internals.h"


// ON WHEN THESE ARE USED

// These functions implement various default priors. The alternative to these
// defaults is using user-specified closures, which only take one parameter
// 'param', and have prior parameters enclosed or hard-coded. If user-specified
// functions are using C/C++ code, we strongly recommend using the native R API
// to distribution functions. For a list of available distributions, see:
// https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions




// ON THE USE OF CLOSURES

// User-specified prior functions are R closures which contain the prior
// parameters, so that the evaluation of the prior is simple and takes a single
// argument 'param'. On the C++ side, using closures would add unwanted
// complexity to the code. Besides, closures would have to be created within the
// movement functions, i.e. every time they are called. Therefore, the C++
// priors used by default are not closures, and take three arguments:

// - param: a Rcpp:List containing parameters

// - config: a Rcpp:List containing parameters for the priors in
// config["prior_xxx"] where 'xxx' is the name of the relevant parameter

// - custom_function: an optional custom prior function (closure); if NULL
// (i.e. R_NilValue in C++) then the basic functions are used


// The prior for the mutation rate 'mu' is an exponential distribution

// [[Rcpp::export(rng = false)]]
double cpp_prior_mu(Rcpp::List param, Rcpp::List config,
		    Rcpp::RObject custom_function = R_NilValue) {

  if (custom_function == R_NilValue) {
    // note: R::dexp is parametrised by scale, not rate
    double scale = 1.0 / Rcpp::as<double>(config["prior_mu"]);

    return R::dexp(Rcpp::as<double>(param["mu"]), scale, true);
 } else {
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(param));
 }

}






// The prior for the reporting probability 'pi' is a beta distribution

// [[Rcpp::export(rng = false)]]
double cpp_prior_pi(Rcpp::List param, Rcpp::List config,
		    Rcpp::RObject custom_function = R_NilValue) {

  if (custom_function == R_NilValue) {
    Rcpp::NumericVector shape = config["prior_pi"];

    return R::dbeta(Rcpp::as<double>(param["pi"]),
			  (double) shape[0],
			  (double) shape[1],
			  true);
  } else {
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(param));
  }

}



// The prior for the trasfer probability 'sigma' is a beta distribution

// [[Rcpp::export(rng = false)]]
double cpp_prior_sigma(Rcpp::List param, Rcpp::List config,
                       Rcpp::RObject custom_function = R_NilValue) {
    
    if (custom_function == R_NilValue) {
        Rcpp::NumericVector shape = config["prior_sigma"];
        
        return R::dbeta(Rcpp::as<double>(param["sigma"]),
                        (double) shape[0],
                        (double) shape[1],
                        true);
    } else {
        Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);
        
        return Rcpp::as<double>(f(param));
    }
    
}



// The prior for the contact reporting coverage 'eps' is a beta distribution

// [[Rcpp::export(rng = false)]]
double cpp_prior_eps(Rcpp::List param, Rcpp::List config,
		    Rcpp::RObject custom_function = R_NilValue) {

  if (custom_function == R_NilValue) {
    Rcpp::NumericVector shape = config["prior_eps"];

    return R::dbeta(Rcpp::as<double>(param["eps"]),
			  (double) shape[0],
			  (double) shape[1],
			  true);
  } else {
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(param));
  }

}






// The prior for the non-transmision contact rate 'lambda' is a beta distribution

// [[Rcpp::export(rng = false)]]
double cpp_prior_lambda(Rcpp::List param, Rcpp::List config,
		    Rcpp::RObject custom_function = R_NilValue) {

  if (custom_function == R_NilValue) {
    Rcpp::NumericVector shape = config["prior_lambda"];

    return R::dbeta(Rcpp::as<double>(param["lambda"]),
			  (double) shape[0],
			  (double) shape[1],
			  true);
  } else {
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);

    return Rcpp::as<double>(f(param));
  }

}






// [[Rcpp::export(rng = false)]]
double cpp_prior_all(Rcpp::List param, Rcpp::List config,
		     Rcpp::RObject custom_functions = R_NilValue
		     ) {
  if (custom_functions == R_NilValue) {

    return cpp_prior_mu(param, config) +
      cpp_prior_pi(param, config) +
      cpp_prior_eps(param, config) +
      cpp_prior_lambda(param, config);

  } else {

    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);

    return cpp_prior_mu(param, config, list_functions["mu"]) +
      cpp_prior_pi(param, config, list_functions["pi"]) +
      cpp_prior_eps(param, config, list_functions["eps"]) +
      cpp_prior_lambda(param, config, list_functions["lambda"]);
  }
}
