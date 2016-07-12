#include <Rcpp.h>
#include <Rmath.h>



/*
 Movement of infection dates are +/- 1 from current states. These movements are currently
 vectorised, i.e. a bunch of dates are proposed all together; this may not be sustainable for
 larger datasets. The non-vectorised option will be slower and speed-up with C/C++ will be more
 substantial then.

 This version differs from the initial R implementation in several points:
 1. all cases are moved
 2. cases are moved one by one
 3. movement for each case is +/- 1 time unit

*/


// [[Rcpp::export("cpp.move.t.inf", rng = true)]]
Rcpp::List cpp_move_t_inf(Rcpp::List data, Rcpp::List param) {
  Rcpp::List new_param(param);
  //Rcpp::IntegerVector t_inf = new_param["t.inf"];
  size_t N = static_cast<size_t>(data["N"]);
  size_t i = 0;

  for (i = 0; i < N; i++) {
    // NEED TO BE ABLE TO CHANGE ELEMENTS OF A LIST HERE
}

  return new_param;
 
}
