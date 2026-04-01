#' Process input pmf for outbreaker_data
#'
#' This function replaces the 0s from a pmf with a minimum probability.
#' 
#' @param dens a vector of probabilities (pmf).
#'
#' @author Cyril Geismar (\email{c.geismar21@@imperial.ac.uk}).
#' 
#' @export
#'
#' @examples
#' 
#' x <- c(0, 0, 0.1, 0.2, 0.3, 0.2, 0, 0.2, 0, 0)
#' sanitize_pmf(x)
#' 

sanitize_pmf <- function(dens){
  
  is_positive <- is.finite(log(dens)) #finds positive values
  
  to_replace <- !is_positive #finds non positive values
  
  min_val <- min(dens[is_positive]) #finds the minimum probability
  
  val_replace <- 1e-4*min_val #generate a smaller probability than the min_val
  
  dens[to_replace] <- val_replace  #replaces non positive values with a tiny probability 
  
  dens <- dens / sum(dens) # rescale density to one
  
  return(dens)
}

