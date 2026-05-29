#' Summary method for bayestsm objects
#'
#' @param object An object of class `"bayestsm"`.
#' @param warmup Number of initial MCMC iterations to discard.
#' @param probs Numeric vector of posterior quantiles.
#' @param ... Additional arguments, currently unused.
#'
#' @method summary bayestsm
#' @export
summary.bayestsm <- function(object, warmup = 0, probs = c(0.025, 0.5, 0.975), ...) {

  if (!inherits(object, "bayestsm")) {
    stop("`object` must be of class 'bayestsm'.", call. = FALSE)
  }

  mod <- object

  make_table <- function(summary_obj, block, mod, probs) {

    quantile_cols <- paste0(probs * 100, "%")

    qtab <- summary_obj$quantiles[, quantile_cols, drop = FALSE]
    qtab <- as.data.frame(qtab, check.names = FALSE)
    qtab$parameter <- rownames(qtab)

    ctab <- mod$convergence$table[
      mod$convergence$table$block == block,
      c("parameter", "R_hat", "ESS"),
      drop = FALSE
    ]

    out <- merge(
      qtab,
      ctab,
      by = "parameter",
      all.x = TRUE,
      sort = FALSE
    )

    rownames(out) <- out$parameter
    out$parameter <- NULL

    out
  }

  summary_X <- summary(mod$par.X)
  summary_S <- summary(mod$par.S)

  table_X <- make_table(summary_X, "X", mod, probs)
  table_S <- make_table(summary_S, "S", mod, probs)

  out <- c(
    "Parameters of x-model",
    capture.output(print(round(table_X, 3))),
    "",
    "Parameters of s-model",
    capture.output(print(round(table_S, 3))),
    "",
    mod$convergence$criteria,
    "",
    sprintf(
      "Total posterior draws saved after thining: %d (total draws: %d)",
      mod$convergence$n_iter * length(mod$par.X),
      mod$thining * mod$convergence$n_iter * length(mod$par.X)
    )
  )

  cat(paste(out, collapse = "\n"))
}

