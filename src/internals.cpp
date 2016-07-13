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
   This function returns the descendents of a given case 'i' in the current ancestries.
*/
// Rcpp::IntegerVector find_descendents(Rcpp::IntegerVector alpha, size_t i) {
//   return match(alpha, i);
// }
