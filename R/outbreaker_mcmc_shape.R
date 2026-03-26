## #' Shape MCMC samples into outputs for outbreaker
## #'
## #' This function shapes MCMC samples for outbreaker into a data.frame with the class \code{outbreaker_chains}.
## #'
## #' @author Thibaut Jombart (\email{thibautjombart@@gmail_com})
## #'
## #' @param param a list of output items as returned by \code{create_param}
## #'
## #' @param data a list of data items as returned by \code{outbreaker_data}
## #'
## #' @export
## #'
outbreaker_mcmc_shape <- function(param, data, config = NULL) {
  ## unfold ancestries ##
  if (!all(vapply(param$alpha, length, integer(1)) == data$N)) {
    stop("some ancestries are missing in the param")
  }
  param$alpha <- matrix(unlist(param$alpha), ncol = data$N, byrow = TRUE)
  colnames(param$alpha) <- paste("alpha", seq_len(data$N), sep = ".")

  ## unfold infection dates ##
  if (!all(vapply(param$t_inf, length, integer(1)) == data$N)) {
    stop("some infection dates are missing in the param")
  }
  param$t_inf <- matrix(unlist(param$t_inf), ncol = data$N, byrow = TRUE)
  colnames(param$t_inf) <- paste("t_inf", seq_len(data$N), sep = ".")

  if (!is.null(data$ctd) || !is.null(data$ctd_timed)) {
    ## unfold epsilon estimates ##
    eps <- matrix(
      unlist(param$eps),
      ncol = length(data$ctd_matrix) + length(data$ctd_timed_matrix),
      byrow = TRUE
    )

    if (ncol(eps) > 1) {
      colnames(eps) <- paste0("eps_", seq_len(ncol(eps)))
    } else {
      colnames(eps) <- 'eps'
    }

    if (!is.null(data$ctd)) {
      ## unfold eta estimates ##
      eta <- matrix(
        unlist(param$eta),
        ncol = length(data$ctd_matrix),
        byrow = TRUE
      )

      if (ncol(eta) > 1) {
        colnames(eta) <- paste0("eta_", seq_len(ncol(eta)))
      } else {
        colnames(eta) <- 'eta'
      }

      ## unfold lambdailon estimates ##
      lambda <- matrix(
        unlist(param$lambda),
        ncol = length(data$ctd_matrix),
        byrow = TRUE
      )

      if (ncol(lambda) > 1) {
        colnames(lambda) <- paste0("lambda_", seq_len(ncol(lambda)))
      } else {
        colnames(lambda) <- 'lambda'
      }
    }
  }

  ## unfold t_onw dates ##
  if (!is.null(data$ctd_timed)) {
    if (!all(vapply(param$t_onw, length, integer(1)) == data$N)) {
      stop("some onward infection dates are missing in the param")
    }
    t_onw <- matrix(unlist(param$t_onw), ncol = data$N, byrow = TRUE)
    t_onw[t_onw == -1000] <- NA
    colnames(t_onw) <- paste("t_onw", seq_len(data$N), sep = ".")

    ## if (!all(vapply(param$place, length, integer(1))==data$N)) {
    ##   stop("some onward infection dates are missing in the param")
    ## }

    ## We don't output the inferred place for now - it would be 3xn additional
    ## columns in the output, mostly of little interest

    ## place <- matrix(unlist(param$place), ncol = data$N, byrow
    ## = TRUE) place[place == 0] <- NA colnames(place) <- paste("place",
    ## seq_len(data$N), sep=".")

    ## unfold tauilon estimates ##
    tau <- matrix(
      unlist(param$tau),
      ncol = length(data$ctd_timed_matrix),
      byrow = TRUE
    )

    colnames(tau) <- paste("tau", seq_len(ncol(tau)), sep = ".")
  }

  ## unfold number of generations ##
  if (!all(vapply(param$kappa, length, integer(1)) == data$N)) {
    stop("some ancestries are missing in the param")
  }
  param$kappa <- matrix(unlist(param$kappa), ncol = data$N, byrow = TRUE)
  colnames(param$kappa) <- paste("kappa", seq_len(data$N), sep = ".")

  ## shape data.frame and convert ##
  param <- data.frame(
    step = param$step,
    post = param$post,
    like = param$like,
    prior = param$prior,
    mu = param$mu,
    pi = param$pi,
    param$alpha,
    param$t_inf,
    param$kappa
  )

  if (!is.null(data$ctd)) {
    param <- cbind(param, eps, lambda, eta)
  }

  if (!is.null(data$ctd_timed)) {
    if (!is.null(data$ctd)) {
      param <- cbind(param, tau, t_onw)
    } else {
      param <- cbind(param, eps, tau, t_onw)
    }
    ##    param <- cbind(param, place)
  }

  names(param) <- gsub("[.]", "_", names(param))

  ## translate integer indexes to case IDs if requested
  if (isTRUE(config$return_ids) && !is.null(data$ids)) {
    ids <- data$ids

    ## rename per-case columns: alpha_1 → alpha_<id>
    alpha_cols <- grep("^alpha_", names(param))
    names(param)[alpha_cols] <- paste0("alpha_", ids)
    t_inf_cols <- grep("^t_inf_", names(param))
    names(param)[t_inf_cols] <- paste0("t_inf_", ids)
    kappa_cols <- grep("^kappa_", names(param))
    names(param)[kappa_cols] <- paste0("kappa_", ids)
    t_onw_cols <- grep("^t_onw_", names(param))
    if (length(t_onw_cols) > 0) {
      names(param)[t_onw_cols] <- paste0("t_onw_", ids)
    }

    ## translate alpha values (the only cross-reference columns)
    for (col in alpha_cols) {
      param[[col]] <- ids[param[[col]]]
    }
  }

  ## output is a data.frame containing all parameters and augmented data, with a dedicated
  ## class (for summary, plotting, etc.)
  attr(param, "ids") <- data$ids
  attr(param, "return_ids") <- isTRUE(config$return_ids)
  attr(param, "min_date") <- data$min_date
  class(param) <- c("outbreaker_chains", "data.frame")
  return(param)
}
