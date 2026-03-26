#' Compute offspring distribution
#'
#' Computes the offspring (secondary case) distribution across posterior
#' samples. For each MCMC step, the reproduction number of every case is
#' computed via \code{\link{get_Ri}}, then pooled into a probability mass
#' function (PMF).
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param burnin The number of iterations to be discarded as burnin.
#' @param raw Logical. If \code{TRUE}, returns a long-format data frame with
#'   columns \code{step}, \code{x}, and \code{y}. If \code{FALSE} (default),
#'   returns summary statistics across steps.
#' @param stats A named list of summary functions applied across steps.
#'   Ignored when \code{raw = TRUE}.
#'
#' @return When \code{raw = FALSE}, a data frame with columns \code{x} (number
#'   of secondary cases), \code{mean}, \code{lwr}, \code{upr}. When
#'   \code{raw = TRUE}, a data frame with columns \code{step}, \code{x},
#'   \code{y}.
#'
#' @author Cyril Geismar (\email{c.geismar21@@imperial.ac.uk}).
#'
#' @seealso \code{\link{get_Ri}} for per-case reproduction numbers;
#'   \code{\link{get_trees}}; \code{\link{outbreaker_chains}}.
#'
#' @importFrom stats aggregate quantile
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
#' get_offspring(out)
#' get_offspring(out, raw = TRUE)
#' }
#'
#' @export
get_offspring <- function(
  out,
  burnin = 0,
  raw = FALSE,
  stats = list(
    mean = mean,
    lwr = function(x) quantile(x, 0.025, na.rm = TRUE),
    upr = function(x) quantile(x, 0.975, na.rm = TRUE)
  )
) {
  ri_long <- get_Ri(out, burnin = burnin, raw = TRUE)
  ri_by_step <- split(ri_long$ri, ri_long$step)

  x_range <- seq(0L, max(ri_long$ri, na.rm = TRUE), by = 1L)
  n_cases <- length(ri_by_step[[1]])
  steps <- unique(ri_long$step)

  off_pmf <- do.call(
    rbind,
    lapply(seq_along(steps), function(j) {
      counts <- tabulate(ri_by_step[[j]] + 1L, nbins = length(x_range))
      data.frame(step = steps[j], x = x_range, y = counts / n_cases)
    })
  )

  if (raw) {
    rownames(off_pmf) <- NULL
    return(off_pmf)
  }

  result <- data.frame(x = x_range)
  for (stat_name in names(stats)) {
    stat_values <- stats::aggregate(y ~ x, off_pmf, stats[[stat_name]])
    result[[stat_name]] <- stat_values$y
  }

  result
}
