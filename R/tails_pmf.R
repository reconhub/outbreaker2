#' Adds tail to pmf containing starting and/or ending 0s
#'
#' This function replaces the tails with 0s with exponential values.
#' 
#' @param dens a vector of probabilities (pmf).
#'
#' @author Cyril Geismar (\email{c.geismar21@@imperial.ac.uk}).
#'
#' @export
#' 
#' 
#' @examples
#' 
#' x <- c(0, 0, 0.1, 0.2, 0.3, 0.2, 0, 0.2, 0, 0)
#' tails_pmf(x)
#' 

tails_pmf <- function(dens){
  
  is_positive <- is.finite(log(dens)) # finds positive values
  
  first_positive <- min(which(is.element(is_positive, TRUE))) # fist positive index
  
  last_positive <- max(which(is.element(is_positive, TRUE))) # last positive index
  
  min_val <- min(dens[is_positive]) #finds minimum value
  
  
  
  ## starting 0s 
  if( isFALSE(is_positive[1]) ){ #if the 1st value is 0
    
    starting_indices <- seq_len(first_positive - 1) # finds starting 0s indices 
    
    val_to_replace <- rev(stats::dexp(starting_indices, 1)) # reverse exponential densities
    
    val_to_replace <- (min_val*1e-4)*(val_to_replace/sum(val_to_replace)) # sum to a value smaller than any other value in dens
    
    dens <- replace(dens, starting_indices, val_to_replace ) # replace satrting 0s with exp densities
    
    warning("starting 0s were replaced with small exponential values")
    
  }
  ## ending 0s
  if( isFALSE(is_positive[length(is_positive)]) ){ #if the last index returns 0
    
    ending_indices <- (last_positive+1) : length(dens) # finds ending 0s indices 
    
    val_to_replace <- stats::dexp(1:length(ending_indices), 1) # reverse exponential densities
    
    val_to_replace <- (min_val*1e-4)*(val_to_replace/sum(val_to_replace)) # sum to a value smaller than any other value in dens
    
    dens <- replace(dens, ending_indices, val_to_replace ) # replace satrting 0s with exp densities
    
    warning("ending 0s were replaced with small exponential values")
    
  }
  
  dens <- dens/sum(dens) #rescale dens to sum to 1
  
  return(dens)
}