#' Compute accuracy of outbreak reconstruction
#'
#' Computes the proportion of correctly assigned infectors in each posterior
#' tree compared to a known true tree. Only useful in simulation studies.
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param true_tree A data frame with \code{from} and \code{to} columns
#'   representing the true transmission tree.
#' @param burnin The number of iterations to be discarded as burnin.
#'
#' @return A numeric vector of accuracy values (0 to 1), one per posterior
#'   sample.
#'
#' @author Cyril Geismar (\email{c.geismar21@@imperial.ac.uk}).
#'
#' @examples
#' \dontrun{
#' data(fake_outbreak)
#' out <- outbreaker(
#'   data = outbreaker_data(
#'     dates = fake_outbreak$onset,
#'     dna = fake_outbreak$dna,
#'     w_dens = fake_outbreak$w
#'   ),
#'   config = create_config(n_iter = 500, sample_every = 50)
#' )
#' true_tree <- data.frame(
#'   from = as.character(fake_outbreak$ances),
#'   to = as.character(seq_along(fake_outbreak$onset))
#' )
#' get_accuracy(out, true_tree)
#' }
#'
#' @export
get_accuracy <- function(out, true_tree, burnin = 0) {
  if (!inherits(out, "outbreaker_chains")) {
    stop("'out' must be an 'outbreaker_chains' object.")
  }
  if (!all(c("from", "to") %in% names(true_tree))) {
    stop("'true_tree' must have 'from' and 'to' columns.")
  }

  trees <- get_trees(out, burnin = burnin)
  true_pairs <- paste(true_tree$from, true_tree$to, sep = "->")

  vapply(
    trees,
    function(df) {
      posterior_pairs <- paste(df$from, df$to, sep = "->")
      sum(posterior_pairs %in% true_pairs) / length(true_pairs)
    },
    numeric(1)
  )
}
