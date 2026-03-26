#' Extract posterior transmission trees
#'
#' Extracts a list of transmission trees from an \code{outbreaker_chains}
#' object, one per posterior sample. Each tree is a data frame with \code{from}
#' and \code{to} columns, optionally augmented with \code{kappa}, \code{t_inf},
#' and user-supplied columns.
#'
#' @param out An object of class \code{outbreaker_chains}.
#' @param burnin The number of iterations to be discarded as burnin.
#' @param kappa Logical. If \code{TRUE}, includes \code{kappa} values.
#' @param t_inf Logical. If \code{TRUE}, includes infection times.
#' @param ... Additional vectors to include as columns. Each vector must have
#'   the same length as the number of cases. The argument name is used as a
#'   column suffix (e.g. \code{date = onset_dates} creates \code{from_date}
#'   and \code{to_date}).
#'
#' @return A list of data frames, one per posterior sample.
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
#' trees <- get_trees(out, date = fake_outbreak$onset)
#' str(trees[[1]])
#' }
#'
#' @export
get_trees <- function(out, burnin = 0, kappa = FALSE, t_inf = FALSE, ...) {
  if (!inherits(out, "outbreaker_chains")) {
    stop("'out' must be an 'outbreaker_chains' object.")
  }

  args <- list(...)
  if (length(args) > 0) {
    stopifnot(all(vapply(args, is.atomic, logical(1))))
  }

  if (burnin > max(out$step)) {
    stop("burnin exceeds the number of steps in 'out'")
  }
  out <- out[out$step > burnin, , drop = FALSE]

  cols <- names(out)
  alpha_cols <- grep("^alpha_", cols, value = TRUE)
  kappa_cols <- grep("^kappa_", cols, value = TRUE)
  t_inf_cols <- grep("^t_inf_", cols, value = TRUE)
  to <- as.character(sub("^alpha_", "", alpha_cols))

  trees <- lapply(seq_len(nrow(out)), function(i) {
    from <- as.character(unlist(out[i, alpha_cols], use.names = FALSE))
    df <- data.frame(from = from, to = to, stringsAsFactors = FALSE)

    if (kappa) {
      kappa_values <- unlist(out[i, kappa_cols], use.names = FALSE)
      names(kappa_values) <- to
      df$from_kappa <- kappa_values[from]
      df$to_kappa <- kappa_values[to]
    }

    if (t_inf) {
      t_inf_values <- unlist(out[i, t_inf_cols], use.names = FALSE)
      names(t_inf_values) <- to
      df$from_t_inf <- t_inf_values[from]
      df$to_t_inf <- t_inf_values[to]
    }

    for (arg in names(args)) {
      vec <- args[[arg]]
      names(vec) <- to
      df[[paste0("from_", arg)]] <- vec[from]
      df[[paste0("to_", arg)]] <- vec[to]
    }

    df
  })

  trees
}
