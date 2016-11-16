
#include <Rcpp.h>



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

//   This function returns a vector of indices of cases which could be infector
//   of 'i' (i.e., their infection dates preceed that of 'i'). Only tricky bit
//   here is keep in mind that 't_inf' is indexed from 0 to N-1, while 'i' and
//   'alpha' (ancestors) are values from 1 to N.

//   Original R code:

// are.possible.alpha <- function(t.inf, i) {
//     if (length(i)>1) {
//         stop("i has a length > 1")
//     }
//     if (any(t.inf[i]==min(t.inf))) {
//         return(NA)
//     }
//     return(which(t.inf < t.inf[i[1]]))
// }

// [[Rcpp::export("cpp.are.possible.ancestors")]]
std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i) {
  size_t n = t_inf.size();
  std::vector<int> out;
  out.reserve(n);
  for (size_t j = 0; j < n; j++) {
    if (t_inf[j] < t_inf[i-1]) { // offset
      out.push_back(j+1); // +1 needed for range 1 ... N
    }
  }
  return out;
}





// ---------------------------

//  This function samples a single value from a vector of integers.

// [[Rcpp::export("cpp.sample1")]]
size_t cpp_sample1(std::vector<int> x) {
  if (x.size() < 1) {
    Rcpp::Rcerr << "Trying to sample from empty vector" << std::endl;
    Rcpp::stop("Trying to sample from empty vector");
  }

  return x[unif_rand() * x.size()];
}




// ---------------------------

//    This function choose a possible infector for case 'i'; 'i' is on the scale
//    1:N

//    Original R version:

// .choose.possible.alpha <- function(t.inf, i) {
//     return(sample(are.possible.alpha(t.inf=t.inf, i=i), 1))
// }

// [[Rcpp::export("cpp.pick.possible.ancestor")]]
size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i) {
  return cpp_sample1(cpp_are_possible_ancestors(t_inf, i));
}






// ---------------------------

// This function returns the descendents of a given case 'i' in the current
// ancestries; 'i' is on the scale 1:N. The output is also on the scale 1:N.

// Original R version:

// find.descendents <- function(param, i) {
//   ## find descendents
//     which(param.current$alpha==i)
//  }

// [[Rcpp::export(cpp.find.descendents)]]
Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i) {
  size_t counter = 0, n = 0;

  // determine size of output vector and create it
  for (size_t j = 0; j < alpha.size(); j++) {
    if (alpha[j] == i) n++;
  }
  
  Rcpp::IntegerVector out(n);

  // fill in output vector
  for (size_t j = 0; j < alpha.size(); j++) {
    if (alpha[j] == i) {
      out[counter++] = j + 1; // offset
    }
  }
  return out;
}







// ---------------------------

// This function returns a vector of indices of cases which are 'local' to a
// case 'i'. Locality is defined as the following set of cases:

// - 'i'
// - the descendents of 'i'
// - 'alpha[i-1]'
// - the descendents of 'alpha[i]' (excluding 'i')

// where 'alpha' is a IntegerVector storing ancestries. Note that 'i' and
// 'alpha' are on the scale 1:N. 

// [[Rcpp::export("cpp.find.local.cases")]]
Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i) {
  // determine descendents of 'i':
  Rcpp::IntegerVector desc_i = cpp_find_descendents(alpha, i);
  size_t n = desc_i.size() + 1; // +1 is to count 'i' itself
  
  // determine descendents of 'alpha[i]':
  Rcpp::IntegerVector desc_alpha_i = cpp_find_descendents(alpha,
							  (size_t) alpha[i-1]);
  if (alpha[i-1] != NA_INTEGER) {
    n += desc_alpha_i.size();
  }

  // create output
  Rcpp::IntegerVector out(n);
  size_t counter = 0;

  // 'i'
  out[counter++] = i;

  // 'descendents of 'i'
  for (size_t j = 0; j < desc_i.size(); j++) {
    out[counter++] = desc_i[j];
  }
  
  if (alpha[i-1] != NA_INTEGER) {
    // alpha[i-1] ...
    out[counter++] = alpha[i-1];
    
    // ... and its descendents
    for (size_t j = 0; j < desc_alpha_i.size(); j++) {
      if ( desc_alpha_i[j] != i) {
	out[counter++] = desc_alpha_i[j];
      }
    }
  }

  return out;
}





// // [[Rcpp::export]]
// std::vector<int> whatever(Rcpp::IntegerVector x) {
//   size_t n = x.size();
//   std::vector<int> foo;
//   foo.reserve(n);
//   for (size_t i; i < n; ++i) {
//     if (condition(x[i])) {
//       foo.push_back(x[i]);
//     }
//   }
//   return foo;
// }

// Rcpp::IntegerVector whatever(Rcpp::IntegerVector x) {
//   size_t n = x.size();
//   std::vector<bool> idx(n);
//   size_t m = 0;
//   for (size_t i = 0; i < n; ++i) {
//     if (condition(x[i])) {
//       ++m;
//       idx[i] = true;
//     }
//   }
//   Rcpp::IntegerVector foo(m);
//   for (size_t i = 0, j = 0; i < n; ++i) {
//     if (idx[i]) {
//       foo[j++] = x[i];
//     }
//   }
//   return foo;
// }









/*
   This function swaps a case 'i' with its infector 'x', so that: 
   - x-> i becomes i->x
   - descendents of i become descendents of x
   - descendents of x become descendents of i

   Original R version:

swap.cases <- function(param, config, i) {
    ## stop if 'i' out of range
    if (i>length(param.current$alpha)) {
        stop("trying to swap ancestry of case ",
             i, " while there are only ",
             length(param.current$alpha), " cases")
    }

    ## find cases for which ancestries can move
    id.ok.to.swap <- which(can.be.swapped(param, config))

    ## find ancestor of 'i'
    x <- param.current$alpha[i]

    ## stop if case 'i' is imported - this should not happen
    if (is.na(x)) {
        warning("trying to swap the ancestry of the imported case ", i)
        return(param)
    }

    ## check that x can be swapped, stop if not
    if (!(x %in% id.ok.to.swap)) {
        return(param)
    }

    ## find indices to swap
    to.be.x <- intersect(which(param.current$alpha==i), id.ok.to.swap)
    to.be.i <- intersect(which(param.current$alpha==x), id.ok.to.swap)

    ## swap 'i' and 'x' in ancestries
    param.current$alpha[to.be.x] <- x
    param.current$alpha[to.be.i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    param.current$alpha[i] <- param.current$alpha[x]

    ## 'i' is now the ancestor of 'x'
    param.current$alpha[x] <- i

    ## swap t.inf
    param.current$t.inf[c(x,i)] <- param.current$t.inf[c(i,x)]

    return(param)
}

*/







// /*
//   This function copies the content if the vector 'a' into the vector 'b'; it throws an error if
//  the sizes of the vectors don't match.
// */

// // [[Rcpp::export]]
// void copy_values(Rcpp::IntegerVector a, Rcpp::IntegerVector b) {
//   size_t N = a.size(), i = 0;

//   if (N != b.size()) {
//     Rcpp::Rcerr << "Trying to copy vectors of different sizes: " << N << " vs " << b.size() << std::endl;
//     Rcpp::stop("Error copying vectors (sizes differ)");
//   }

//   for (i = 0; i <  N; i++) {
//     b[i] = a[i];
//   }
// }
