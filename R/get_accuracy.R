#' Compute accuracy of outbreak reconstruction
#'
#' Compares posterior transmission trees to a known true ancestry vector.
#' Correctly inferred imports (both true and inferred infectors are \code{NA})
#' count as correct.
#'
#' Two modes are available:
#' \itemize{
#'   \item \code{by_case = TRUE} (default): for each case, what proportion of
#'     MCMC steps correctly identified its infector? This tells you which cases
#'     are well-resolved.
#'   \item \code{by_case = FALSE}: for each MCMC step, what proportion of cases
#'     had their infector correctly assigned? This tells you how good the
#'     overall tree is across the posterior.
#' }
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param ances A vector of true infectors, where \code{ances[i]} is the
#'   infector of case \code{i} and \code{NA} indicates an imported case. Must
#'   be the same length and order as the cases used in \code{\link{outbreaker}}.
#' @param burnin The number of iterations to be discarded as burnin.
#' @param by_case Logical. If \code{TRUE} (default), returns one accuracy
#'   value per case. If \code{FALSE}, returns one accuracy value per MCMC step.
#'
#' @return A numeric vector of values between 0 and 1. Named by case
#'   identifier when \code{by_case = TRUE}.
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
#' get_accuracy(out, fake_outbreak$ances)
#' get_accuracy(out, fake_outbreak$ances, by_case = FALSE)
#' }
#'
#' @export
get_accuracy <- function(out, ances, burnin = 0, by_case = TRUE) {
  if (!inherits(out, "outbreaker_chains")) {
    stop("'out' must be an 'outbreaker_chains' object.")
  }

  if (burnin > max(out$step)) {
    stop("burnin exceeds the number of steps in 'out'")
  }
  out <- out[out$step > burnin, , drop = FALSE]

  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  case_ids <- sub("^alpha_", "", alpha_cols)
  alpha_mat <- as.matrix(out[, alpha_cols])

  ances <- as.character(ances)

  correct <- t(apply(alpha_mat, 1, function(row) {
    mapply(identical, as.character(row), ances)
  }))

  agg <- if (by_case) colMeans else rowMeans
  acc <- agg(correct, na.rm = TRUE)
  if (by_case) {
    names(acc) <- case_ids
  }

  return(acc)
}
