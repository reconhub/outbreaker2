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



// // [[Rcpp::export]]
// void add_one_to_x(Rcpp::List a) {
// // Rcpp::List add_one_to_x(Rcpp::List a) {
  
//   //  Rcpp::NumericVector new_x = clone(a["x"]);
//   Rcpp::IntegerVector x = a["x"];
//   Rcpp::IntegerVector new_x = x;
//   //  Rcpp::List out = clone(a);
  
//   size_t i, N = new_x.size();
//   for(i = 0; i < N; i++) {
//     new_x[i] += 1;
//   }
  
//   //x = new_x; // this does not work
//   a["x"] = new_x;
//   // return out;
// }

/* 
   This function returns the descendents of a given case 'i' in the current ancestries.
*/
// Rcpp::IntegerVector find_descendents(Rcpp::IntegerVector alpha, size_t i) {
//   return match(alpha, i);
// }
