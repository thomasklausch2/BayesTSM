#' handle_bayestsm_convergence
#' @noRd
handle_bayestsm_convergence <- function(
    ret,
    update_till_convergence,
    chains,
    min_R,
    min_eff,
    maxit,
    burnin,
    mc_update,
    silent,
    fix.sigma.X = FALSE,
    fix.sigma.S = FALSE
) {
  make_return <- function(ret, needs_update, mc_update_next) {
    list(
      ret = ret,
      convergence = if (!is.null(ret$convergence)) ret$convergence else NULL,
      needs_update = needs_update,
      mc_update_next = mc_update_next
    )
  }

  if (!isTRUE(update_till_convergence)) {
    return(make_return(
      ret = ret,
      needs_update = FALSE,
      mc_update_next = NULL
    ))
  }

  silent <- isTRUE(silent)

  warn_short <- function(msg, silent) {
    if (!silent) {
      warning(msg, call. = FALSE)
    }
  }

  message_short <- function(msg, silent) {
    if (!silent) {
      message(msg)
    }
  }

  as_count <- function(x, name, allow_null = FALSE, min_value = 1L) {
    if (allow_null && is.null(x)) {
      return(NULL)
    }

    if (
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x)
    ) {
      stop(sprintf("%s must be a single finite number.", name), call. = FALSE)
    }

    x <- as.integer(round(x))

    if (is.na(x) || x < min_value) {
      stop(
        sprintf("%s must be >= %d after rounding.", name, min_value),
        call. = FALSE
      )
    }

    x
  }

  as_positive_number <- function(x, name) {
    if (
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0
    ) {
      stop(sprintf("%s must be a single positive number.", name), call. = FALSE)
    }

    x
  }

  as_logical1 <- function(x, name) {
    if (
      !is.logical(x) ||
      length(x) != 1L ||
      is.na(x)
    ) {
      stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
    }

    x
  }

  if (is.null(ret$par.X) || is.null(ret$par.S)) {
    stop("ret must contain par.X and par.S.", call. = FALSE)
  }

  chains <- as_count(
    x = chains,
    name = "chains",
    allow_null = FALSE,
    min_value = 1L
  )

  if (chains < 2L) {
    stop(
      "update_till_convergence = TRUE requires at least 2 chains.",
      call. = FALSE
    )
  }

  if (length(ret$par.X) != chains || length(ret$par.S) != chains) {
    stop(
      "chains must match length(ret$par.X) and length(ret$par.S).",
      call. = FALSE
    )
  }

  maxit <- as_count(
    x = maxit,
    name = "maxit",
    allow_null = TRUE,
    min_value = 1L
  )

  mc_update <- as_count(
    x = mc_update,
    name = "mc_update",
    allow_null = FALSE,
    min_value = 1L
  )

  burnin <- as_count(
    x = burnin,
    name = "burnin",
    allow_null = FALSE,
    min_value = 0L
  )

  min_R <- as_positive_number(min_R, "min_R")
  min_eff <- as_positive_number(min_eff, "min_eff")

  fix.sigma.X <- as_logical1(fix.sigma.X, "fix.sigma.X")
  fix.sigma.S <- as_logical1(fix.sigma.S, "fix.sigma.S")

  get_par_names <- function(x, mats) {
    pnames <- tryCatch(
      coda::varnames(x),
      error = function(e) NULL
    )

    if (is.null(pnames) || length(pnames) == 0L) {
      pnames <- colnames(mats[[1L]])
    }

    if (is.null(pnames)) {
      return(character(0L))
    }

    pnames <- as.character(pnames)
    pnames[!is.na(pnames) & nzchar(pnames)]
  }

  make_block <- function(x, block_name, burnin, silent) {
    mats <- lapply(x, as.matrix)

    n_iter <- vapply(mats, nrow, integer(1L))

    if (length(unique(n_iter)) != 1L) {
      stop(
        sprintf("All chains in block %s must have the same length.", block_name),
        call. = FALSE
      )
    }

    n_iter <- n_iter[[1L]]

    rows <- if (burnin < n_iter) {
      seq.int(from = burnin + 1L, to = n_iter)
    } else {
      integer(0L)
    }

    if (length(rows) < 2L) {
      warn_short(
        sprintf("Too few post-burnin draws for block %s.", block_name),
        silent = silent
      )
    }

    list(
      mats = mats,
      pnames = get_par_names(x = x, mats = mats),
      rows = rows,
      n_iter = n_iter
    )
  }

  make_draw_matrix <- function(mats, rows, param_name) {
    vapply(
      X = mats,
      FUN = function(chain_mat) {
        as.numeric(chain_mat[rows, param_name, drop = TRUE])
      },
      FUN.VALUE = numeric(length(rows))
    )
  }

  drop_fixed_sigma <- function(block, block_name, fix_sigma, silent) {
    block$fixed_pnames <- character(0L)

    if (!isTRUE(fix_sigma)) {
      return(block)
    }

    if (length(block$pnames) == 0L) {
      warn_short(
        sprintf(
          "Could not exclude fixed sigma in block %s because no parameter names were found.",
          block_name
        ),
        silent = silent
      )
      return(block)
    }

    sigma_name <- block$pnames[length(block$pnames)]

    if (length(block$rows) > 0L) {
      sigma_draws <- tryCatch(
        make_draw_matrix(
          mats = block$mats,
          rows = block$rows,
          param_name = sigma_name
        ),
        error = function(e) NULL
      )

      if (!is.null(sigma_draws) && length(sigma_draws) > 0L) {
        tol <- sqrt(.Machine$double.eps)
        is_constant <- all(abs(sigma_draws - sigma_draws[1L, 1L]) <= tol)

        if (!is_constant) {
          warn_short(
            sprintf(
              "Parameter %s in block %s was marked as fixed but is not constant across post-burnin draws.",
              sigma_name,
              block_name
            ),
            silent = silent
          )
        }
      }
    }

    block$pnames <- block$pnames[-length(block$pnames)]
    block$fixed_pnames <- sigma_name

    block
  }

  safe_diag <- function(block, block_name, fun, fun_label, silent) {
    out <- stats::setNames(
      rep(NA_real_, length(block$pnames)),
      block$pnames
    )

    if (length(block$rows) < 2L || length(block$pnames) == 0L) {
      return(out)
    }

    for (param_name in block$pnames) {
      has_param <- vapply(
        X = block$mats,
        FUN = function(chain_mat) param_name %in% colnames(chain_mat),
        FUN.VALUE = logical(1L)
      )

      if (!all(has_param)) {
        warn_short(
          sprintf("Missing parameter %s in block %s.", param_name, block_name),
          silent = silent
        )
        next
      }

      x_p <- tryCatch(
        make_draw_matrix(
          mats = block$mats,
          rows = block$rows,
          param_name = param_name
        ),
        error = function(e) {
          warn_short(
            sprintf("Could not build draws for %s:%s.", block_name, param_name),
            silent = silent
          )
          NULL
        }
      )

      if (is.null(x_p) || nrow(x_p) < 2L || ncol(x_p) < 2L) {
        next
      }

      out[param_name] <- tryCatch(
        as.numeric(fun(x_p)),
        error = function(e) {
          warn_short(
            sprintf("%s failed for %s:%s.", fun_label, block_name, param_name),
            silent = silent
          )
          NA_real_
        }
      )
    }

    out
  }

  meets_upper <- function(x, limit) {
    length(x) > 0L && all(is.finite(x)) && all(x <= limit)
  }

  meets_lower <- function(x, limit) {
    length(x) > 0L && all(is.finite(x)) && all(x >= limit)
  }

  get_thinning <- function(ret) {
    thinning <- if (!is.null(ret$thining)) {
      ret$thining
    } else {
      ret$thinning
    }

    if (
      is.null(thinning) ||
      !is.numeric(thinning) ||
      length(thinning) != 1L ||
      is.na(thinning) ||
      !is.finite(thinning) ||
      thinning <= 0
    ) {
      return(1L)
    }

    as.integer(round(thinning))
  }

  block_x <- make_block(
    x = ret$par.X,
    block_name = "X",
    burnin = burnin,
    silent = silent
  )

  block_s <- make_block(
    x = ret$par.S,
    block_name = "S",
    burnin = burnin,
    silent = silent
  )

  if (block_x$n_iter != block_s$n_iter) {
    stop(
      "Blocks X and S must have the same number of stored iterations.",
      call. = FALSE
    )
  }

  block_x <- drop_fixed_sigma(
    block = block_x,
    block_name = "X",
    fix_sigma = fix.sigma.X,
    silent = silent
  )

  block_s <- drop_fixed_sigma(
    block = block_s,
    block_name = "S",
    fix_sigma = fix.sigma.S,
    silent = silent
  )

  g_x <- safe_diag(
    block = block_x,
    block_name = "X",
    fun = posterior::rhat,
    fun_label = "R-hat",
    silent = silent
  )

  g_s <- safe_diag(
    block = block_s,
    block_name = "S",
    fun = posterior::rhat,
    fun_label = "R-hat",
    silent = silent
  )

  eff_x <- safe_diag(
    block = block_x,
    block_name = "X",
    fun = posterior::ess_mean,
    fun_label = "ESS",
    silent = silent
  )

  eff_s <- safe_diag(
    block = block_s,
    block_name = "S",
    fun = posterior::ess_mean,
    fun_label = "ESS",
    silent = silent
  )

  converged <-
    meets_upper(g_x, min_R) &&
    meets_upper(g_s, min_R) &&
    meets_lower(eff_x, min_eff) &&
    meets_lower(eff_s, min_eff)

  conv_tab <- format_convergence_table(
    g_x = g_x,
    g_s = g_s,
    eff_x = eff_x,
    eff_s = eff_s
  )

  fixed_pnames <- list(
    X = block_x$fixed_pnames,
    S = block_s$fixed_pnames
  )

  n_fixed <- length(unlist(fixed_pnames, use.names = FALSE))

  if (n_fixed > 0L) {
    criteria_text <- sprintf(
      paste0(
        "Convergence criteria: R-hat <= %.3f and ESS >= %.1f ",
        "for sampled parameters; fixed sigma parameters excluded."
      ),
      min_R,
      min_eff
    )
  } else {
    criteria_text <- sprintf(
      "Convergence criteria: R-hat <= %.3f and ESS >= %.1f",
      min_R,
      min_eff
    )
  }

  n_iter <- block_x$n_iter

  ret$convergence <- list(
    g_x = g_x,
    g_s = g_s,
    eff_x = eff_x,
    eff_s = eff_s,
    table = conv_tab,
    min_R = min_R,
    min_eff = min_eff,
    burnin = burnin,
    criteria = criteria_text,
    n_iter = n_iter,
    n_chains = chains,
    fixed_pnames = fixed_pnames,
    converged = converged
  )

  tab_txt <- paste(
    capture.output(print(conv_tab, row.names = FALSE)),
    collapse = "\n"
  )

  if (converged) {
    message_short(
      sprintf(
        paste0(
          "Converged after %d stored iter/chain.\n",
          "%s\n",
          "%s"
        ),
        n_iter,
        criteria_text,
        tab_txt
      ),
      silent = silent
    )

    return(make_return(
      ret = ret,
      needs_update = FALSE,
      mc_update_next = NULL
    ))
  }

  if (!is.null(maxit) && n_iter >= maxit) {
    message_short(
      sprintf(
        paste0(
          "Stopped at maxit = %d without convergence.\n",
          "%s\n",
          "%s"
        ),
        maxit,
        criteria_text,
        tab_txt
      ),
      silent = silent
    )

    return(make_return(
      ret = ret,
      needs_update = FALSE,
      mc_update_next = NULL
    ))
  }

  if (is.null(maxit)) {
    mc_update_next <- mc_update

    message_short(
      sprintf(
        paste0(
          "Not converged after %d stored iter/chain; update +%d raw.\n",
          "%s\n",
          "%s\n\n"
        ),
        n_iter,
        mc_update_next,
        criteria_text,
        tab_txt
      ),
      silent = silent
    )
  } else {
    thinning <- get_thinning(ret = ret)
    remaining_stored_iter <- maxit - n_iter
    remaining_raw_iter <- remaining_stored_iter * thinning
    mc_update_next <- min(mc_update, remaining_raw_iter)

    message_short(
      sprintf(
        paste0(
          "Not converged after %d/%d stored iter/chain; update +%d raw.\n",
          "%s\n",
          "%s\n\n"
        ),
        n_iter,
        maxit,
        mc_update_next,
        criteria_text,
        tab_txt
      ),
      silent = silent
    )
  }

  make_return(
    ret = ret,
    needs_update = TRUE,
    mc_update_next = mc_update_next
  )
}
