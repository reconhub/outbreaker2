#' Compute serial interval distribution
#'
#' Computes the serial interval distribution across posterior samples. The
#' serial interval is the time between onset dates in each infector-infectee
#' pair. For each MCMC step, transmission pairs are identified from the
#' \code{alpha} columns and onset-to-onset delays are computed using the dates
#' in the \code{outbreaker_data} object.
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param data An \code{outbreaker_data} object containing the onset dates.
#' @param burnin The number of iterations to be discarded as burnin.
#' @param raw Logical. If \code{TRUE}, returns a long-format data frame with
#'   columns \code{step}, \code{x}, and \code{y}. If \code{FALSE} (default),
#'   returns summary statistics across steps.
#' @param stats A named list of summary functions applied across steps.
#'   Ignored when \code{raw = TRUE}.
#'
#' @return When \code{raw = FALSE}, a data frame with columns \code{x} (serial
#'   interval value), \code{mean}, \code{lwr}, \code{upr}. When
#'   \code{raw = TRUE}, a data frame with columns \code{step}, \code{x},
#'   \code{y}.
#'
#' @author Cyril Geismar (\email{c.geismar21@@imperial.ac.uk}).
#'
#' @seealso \code{\link{get_trees}}, \code{\link{outbreaker_chains}}.
#'
#' @importFrom stats aggregate quantile
#'
#' @examples
#' \dontrun{
#' data(fake_outbreak)
#' dat <- outbreaker_data(
#'   dates = fake_outbreak$onset,
#'   dna = fake_outbreak$dna,
#'   w_dens = fake_outbreak$w
#' )
#' out <- outbreaker(data = dat,
#'   config = create_config(n_iter = 500, sample_every = 50)
#' )
#' get_si(out, dat)
#' get_si(out, dat, raw = TRUE)
#' }
#'
#' @export
get_si <- function(
  out,
  data,
  burnin = 0,
  raw = FALSE,
  stats = list(
    mean = mean,
    lwr = function(x) quantile(x, 0.025, na.rm = TRUE),
    upr = function(x) quantile(x, 0.975, na.rm = TRUE)
  )
) {

  if (!inherits(out, "outbreaker_chains")) {
    stop("'out' must be an 'outbreaker_chains' object.")
  }
  if (is.null(data$dates)) {
    stop("'data' must contain onset dates ('dates').")
  }

  if (burnin > max(out$step)) {
    stop("burnin exceeds the number of steps in 'out'")
  }
  out <- out[out$step > burnin, , drop = FALSE]

  dates <- data$dates
  names(dates) <- data$ids

  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  case_ids <- sub("^alpha_", "", alpha_cols)

  si_list <- lapply(seq_len(nrow(out)), function(i) {
    alpha <- unlist(out[i, alpha_cols], use.names = FALSE)
    has_ancestor <- !is.na(alpha)
    from_ids <- as.character(alpha[has_ancestor])
    to_ids <- case_ids[has_ancestor]
    as.numeric(dates[to_ids] - dates[from_ids])
  })

  all_si <- unlist(si_list)
  x_range <- seq(floor(min(all_si, na.rm = TRUE)),
                 ceiling(max(all_si, na.rm = TRUE)),
                 by = 1)
  steps <- out$step

  si_pmf <- do.call(
    rbind,
    lapply(seq_along(si_list), function(j) {
      si <- si_list[[j]]
      counts <- tabulate(match(si, x_range), nbins = length(x_range))
      data.frame(step = steps[j], x = x_range, y = counts / length(si))
    })
  )

  if (raw) {
    rownames(si_pmf) <- NULL
    return(si_pmf)
  }

  result <- data.frame(x = x_range)
  for (stat_name in names(stats)) {
    stat_values <- stats::aggregate(y ~ x, si_pmf, stats[[stat_name]])
    result[[stat_name]] <- stat_values$y
  }

  result
}
