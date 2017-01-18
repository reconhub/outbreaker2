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

// User-specified prior functions are closures which contain the prior
// parameters, so that the evaluation of the prior is simple and takes a single
// argument 'param'. On the C++ side, using closures would add unwanted
// complexity to the code. Besides, closures would have to be created within the
// movement functions, i.e. every time they are called. Therefore, the C++
// priors used by default are not closures, and take two arguments:

// - param: a Rcpp:List containing parameters

// - config: a Rcpp:List containing parameters for the priors in
// config["prior_xxx"] where 'xxx' is the name of the relevant parameter



// The prior for the mutation rate 'mu' is an exponential distribution

double cpp_prior_mu(Rcpp::List param, Rcpp::List config) {
  // note: R::dexp is parametrised by scale, not rate

  double scale = 1.0 / Rcpp::as<double>(config["prior_mu"]);

  double out = R::dexp(Rcpp::as<double>(param["mu"]), scale, true);

  return out;
  
}





// The prior for the reporting probability 'pi' is an beta distribution

double cpp_prior_pi(Rcpp::List param, Rcpp::List config) {

  Rcpp::NumericVector shape = config["prior_pi"];

  double out = R::dbeta(Rcpp::as<double>(param["pi"]), 
			(double) shape[0], 
			(double) shape[1], 
			true);
  
  return out;
}



