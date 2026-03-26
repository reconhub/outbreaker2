#' Compute entropy of infector assignments
#'
#' Computes the Shannon entropy of inferred infectors for each case across
#' posterior samples. Entropy quantifies uncertainty in infector assignment:
#' 0 means complete certainty, 1 means maximum uncertainty.
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param burnin The number of iterations to be discarded as burnin.
#' @param normalise Logical. If \code{TRUE} (default), entropy is normalised
#'   to [0, 1] by dividing by \code{log(K)} where \code{K} is the number of
#'   distinct inferred infectors.
#'
#' @return A named numeric vector of entropy values, one per case.
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
#' get_entropy(out)
#' }
#'
#' @export
get_entropy <- function(out, burnin = 0, normalise = TRUE) {
  if (!inherits(out, "outbreaker_chains")) {
    stop("'out' must be an 'outbreaker_chains' object.")
  }

  if (burnin > max(out$step)) {
    stop("burnin exceeds the number of steps in 'out'")
  }
  out <- out[out$step > burnin, , drop = FALSE]

  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  ent <- vapply(
    out[alpha_cols],
    .entropy,
    normalise = normalise,
    FUN.VALUE = numeric(1)
  )
  names(ent) <- sub("^alpha_", "", alpha_cols)
  ent
}


## internal: Shannon entropy of a categorical vector
.entropy <- function(x, normalise = TRUE) {
  p <- table(as.character(x))
  p <- p / sum(p)
  H <- -sum(p * log(p))

  if (!normalise) {
    return(H)
  }

  K <- length(p)
  if (K <= 1) {
    return(0)
  }
  H / log(K)
}
