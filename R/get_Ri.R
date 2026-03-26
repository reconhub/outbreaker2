#' Compute case reproduction numbers
#'
#' Computes the number of secondary infections caused by each case across
#' posterior samples. For each MCMC step, counts how many times each case
#' appears as an infector in the \code{alpha} columns.
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param burnin The number of iterations to be discarded as burnin.
#' @param raw Logical. If \code{TRUE}, returns a long-format data frame with
#'   columns \code{step}, \code{case}, and \code{ri}. If \code{FALSE}
#'   (default), returns summary statistics per case.
#' @param stats A named list of summary functions applied across steps.
#'   Ignored when \code{raw = TRUE}.
#'
#' @return When \code{raw = FALSE}, a data frame with columns \code{case},
#'   \code{mean}, \code{lwr}, \code{upr}. When \code{raw = TRUE}, a data
#'   frame with columns \code{step}, \code{case}, \code{ri}.
#'
#' @author Cyril Geismar (\email{c.geismar21@@imperial.ac.uk}).
#'
#' @seealso \code{\link{get_offspring}}, \code{\link{get_trees}},
#'   \code{\link{outbreaker_chains}}.
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
#' get_Ri(out)
#' get_Ri(out, raw = TRUE)
#' }
#'
#' @export
get_Ri <- function(
  out,
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

  if (burnin > max(out$step)) {
    stop("burnin exceeds the number of steps in 'out'")
  }
  out <- out[out$step > burnin, , drop = FALSE]

  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  case_ids <- sub("^alpha_", "", alpha_cols)
  alpha_mat <- as.matrix(out[, alpha_cols])

  Ri_mat <- t(apply(alpha_mat, 1, function(row) {
    as.integer(table(factor(row, levels = case_ids)))
  }))
  colnames(Ri_mat) <- case_ids

  ri_long <- data.frame(
    step = rep(out$step, each = length(case_ids)),
    case = rep(case_ids, times = nrow(out)),
    ri = as.vector(t(Ri_mat)),
    stringsAsFactors = FALSE
  )

  if (raw) {
    rownames(ri_long) <- NULL
    return(ri_long)
  }

  result <- data.frame(case = case_ids, stringsAsFactors = FALSE)
  for (stat_name in names(stats)) {
    stat_values <- stats::aggregate(ri ~ case, ri_long, stats[[stat_name]])
    result[[stat_name]] <- stat_values$ri[match(case_ids, stat_values$case)]
  }

  result
}
