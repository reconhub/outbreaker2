#include <Rcpp.h>

/*
  This function copies the content if the vector 'a' into the vector 'b'; it throws an error if
 the sizes of the vectors don't match.
*/

// [[Rcpp::export]]
void copy_values(Rcpp::IntegerVector a, Rcpp::IntegerVector b) {
  size_t N = a.size(), i = 0;

  if (N != b.size()) {
    Rcpp::Rcerr << "Trying to copy vectors of different sizes: " << N << " vs " << b.size() << std::endl;
    Rcpp::stop("Error copying vectors (sizes differ)");
  }

  for (i = 0; i <  N; i++) {
    b[i] = a[i];
  }
}







/*
  This function returns a vector of indices of cases which could be infector of 'i' (i.e., their
  infection dates preceed that of 'i').

  Original R code:

are.possible.alpha <- function(t.inf, i) {
    if (length(i)>1) {
        stop("i has a length > 1")
    }
    if (any(t.inf[i]==min(t.inf))) {
        return(NA)
    }
    return(which(t.inf < t.inf[i[1]]))
}

*/







/*
   This function choose a possible infector for case 'i'.

   Original R version:

choose.possible.alpha <- function(t.inf, i) {
    return(sample(are.possible.alpha(t.inf=t.inf, i=i), 1))
}

*/


/*
   This function returns the descendents of a given case 'i' in the current ancestries.

   Original R version:

find.descendents <- function(param, i) {
    ## find descendents
    which(param$current.alpha==i)
}


*/
// Rcpp::IntegerVector find_descendents(Rcpp::IntegerVector alpha, size_t i) {
//   return match(alpha, i);
// }








/*
   This function swaps a case 'i' with its infector 'x', so that: 
   - x-> i becomes i->x
   - descendents of i become descendents of x
   - descendents of x become descendents of i

   Original R version:

swap.cases <- function(param, config, i) {
    ## stop if 'i' out of range
    if (i>length(param$current.alpha)) {
        stop("trying to swap ancestry of case ",
             i, " while there are only ",
             length(param$current.alpha), " cases")
    }

    ## find cases for which ancestries can move
    id.ok.to.swap <- which(can.be.swapped(param, config))

    ## find ancestor of 'i'
    x <- param$current.alpha[i]

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
    to.be.x <- intersect(which(param$current.alpha==i), id.ok.to.swap)
    to.be.i <- intersect(which(param$current.alpha==x), id.ok.to.swap)

    ## swap 'i' and 'x' in ancestries
    param$current.alpha[to.be.x] <- x
    param$current.alpha[to.be.i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    param$current.alpha[i] <- param$current.alpha[x]

    ## 'i' is now the ancestor of 'x'
    param$current.alpha[x] <- i

    ## swap t.inf
    param$current.t.inf[c(x,i)] <- param$current.t.inf[c(i,x)]

    return(param)
}

*/
