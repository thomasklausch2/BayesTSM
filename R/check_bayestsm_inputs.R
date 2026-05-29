#' Check and standardize inputs for bayestsm()
#'
#' Validates and standardizes the observed data, covariate matrices, weights,
#' and proposal specification used by bayestsm(). Incomplete observations are
#' removed with a warning.
#'
#' @param d Numeric/integer vector with values in {1, 2, 3}.
#' @param L Numeric vector of left interval bounds.
#' @param R Numeric vector of right interval bounds. May contain Inf for
#'   right-censored observations.
#' @param Z.X Optional covariate matrix/data.frame for the X model.
#' @param Z.S Optional covariate matrix/data.frame for the S model. If NULL,
#'   defaults to Z.X.
#' @param w Optional weight vector. Must have length 1 or length(d).
#' @param prop.sd Proposal specification. Must be supplied unless prev.run is
#'   not NULL. May be a positive scalar or a numeric square matrix.
#' @param prev.run Optional previous bayestsm fit.
#'
#' @return A list with cleaned and standardized inputs:
#'   \itemize{
#'     \item d, L, R
#'     \item Z.X, Z.S
#'     \item Z1.X, Z1.S
#'     \item p.X, p.S, p1.X, p1.S
#'     \item w
#'     \item prop.sd
#'     \item omitted
#'   }
#'
#' @noRd
check_bayestsm_inputs <- function(d,
                                  L,
                                  R,
                                  Z.X = NULL,
                                  Z.S = NULL,
                                  w = 1,
                                  prop.sd = NULL,
                                  prev.run = NULL) {

  if (is.null(d) || is.null(L) || is.null(R)) {
    stop("Arguments 'd', 'L', and 'R' must all be supplied.", call. = FALSE)
  }

  # if (is.null(prop.sd) && is.null(prev.run)) {
  #   stop("Argument 'prop.sd' must be supplied unless 'prev.run' is provided.",
  #        call. = FALSE)
  # }

  d <- as.vector(d)
  L <- as.vector(L)
  R <- as.vector(R)

  n <- length(d)

  if (length(L) != n || length(R) != n) {
    stop("Arguments 'd', 'L', and 'R' must have the same length.",
         call. = FALSE)
  }

  if (n == 0L) {
    stop("No observations supplied.", call. = FALSE)
  }

  if (is.null(Z.S)) {
    Z.S <- Z.X
  }

  if (!is.null(Z.X)) {
    Z.X <- as.matrix(Z.X)
    if (nrow(Z.X) != n) {
      stop("Argument 'Z.X' must have the same number of rows as length(d).",
           call. = FALSE)
    }
    if (!is.numeric(Z.X)) {
      stop("Argument 'Z.X' must be numeric or coercible to a numeric matrix.",
           call. = FALSE)
    }
    storage.mode(Z.X) <- "double"
    if (is.null(colnames(Z.X))) {
      colnames(Z.X) <- paste0("ZX.", seq_len(ncol(Z.X)))
    }
  }

  if (!is.null(Z.S)) {
    Z.S <- as.matrix(Z.S)
    if (nrow(Z.S) != n) {
      stop("Argument 'Z.S' must have the same number of rows as length(d).",
           call. = FALSE)
    }
    if (!is.numeric(Z.S)) {
      stop("Argument 'Z.S' must be numeric or coercible to a numeric matrix.",
           call. = FALSE)
    }
    storage.mode(Z.S) <- "double"
    if (is.null(colnames(Z.S))) {
      colnames(Z.S) <- paste0("ZS.", seq_len(ncol(Z.S)))
    }
  }

  if (length(w) == 1L) {
    w <- rep(w, n)
  } else {
    w <- as.vector(w)
  }

  if (length(w) != n) {
    stop("Argument 'w' must have length 1 or length(d).", call. = FALSE)
  }

  cc <- complete.cases(d, L, R, w) &
    (if (is.null(Z.X)) rep(TRUE, n) else complete.cases(Z.X)) &
    (if (is.null(Z.S)) rep(TRUE, n) else complete.cases(Z.S))

  omitted <- which(!cc)

  if (length(omitted) > 0L) {
    warning(
      sprintf(
        "Removed %d observation(s) with missing values in d, L, R, w, Z.X, or Z.S.",
        length(omitted)
      ),
      call. = FALSE
    )
    d <- d[cc]
    L <- L[cc]
    R <- R[cc]
    w <- w[cc]
    if (!is.null(Z.X)) {
      Z.X <- Z.X[cc, , drop = FALSE]
    }
    if (!is.null(Z.S)) {
      Z.S <- Z.S[cc, , drop = FALSE]
    }
  }

  n <- length(d)

  if (n == 0L) {
    stop("No complete observations remain after removing missing values.",
         call. = FALSE)
  }

  if (!is.numeric(d)) {
    stop("Argument 'd' must be numeric/integer with values 1, 2, or 3.",
         call. = FALSE)
  }

  if (!is.numeric(L) || !is.numeric(R)) {
    stop("Arguments 'L' and 'R' must be numeric.", call. = FALSE)
  }

  if (!is.numeric(w)) {
    stop("Argument 'w' must be numeric.", call. = FALSE)
  }

  if (any(!d %in% c(1, 2, 3))) {
    stop("Argument 'd' may only contain the values 1, 2, and 3.",
         call. = FALSE)
  }

  if (any(abs(d - round(d)) > .Machine$double.eps^0.5)) {
    stop("Argument 'd' must contain integer-coded values 1, 2, and 3.",
         call. = FALSE)
  }

  d <- as.integer(d)

  if (any(!is.finite(L))) {
    stop("All values of 'L' must be finite.", call. = FALSE)
  }

  if (any(L < 0)) {
    stop("All values of 'L' must be >= 0.", call. = FALSE)
  }

  if (any(!is.finite(R) & R != Inf)) {
    stop("Argument 'R' may only be finite or equal to Inf.", call. = FALSE)
  }

  if (any(R <= L, na.rm = TRUE)) {
    stop("All observations must satisfy R > L.", call. = FALSE)
  }

  if (any(R == Inf & d != 1L)) {
    stop("If R = Inf, then d must be 1.", call. = FALSE)
  }

  if (any(is.finite(R) & d == 1L)) {
    stop("If R is finite, then d must be 2 or 3.", call. = FALSE)
  }

  if (any(!is.finite(w))) {
    stop("All weights in 'w' must be finite.", call. = FALSE)
  }

  if (any(w <= 0)) {
    stop("All weights in 'w' must be > 0.", call. = FALSE)
  }

  if (!is.null(Z.X) && any(!is.finite(Z.X))) {
    stop("All values in 'Z.X' must be finite.", call. = FALSE)
  }

  if (!is.null(Z.S) && any(!is.finite(Z.S))) {
    stop("All values in 'Z.S' must be finite.", call. = FALSE)
  }

  if (!is.null(Z.X) && ncol(Z.X) == 0L) {
    Z.X <- NULL
  }

  if (!is.null(Z.S) && ncol(Z.S) == 0L) {
    Z.S <- NULL
  }

  if (is.null(Z.X)) {
    p.X <- 0L
    Z1.X <- matrix(1, nrow = n, ncol = 1L)
    colnames(Z1.X) <- "Intercept"
  } else {
    p.X <- ncol(Z.X)
    Z1.X <- cbind(Intercept = 1, Z.X)
  }

  if (is.null(Z.S)) {
    p.S <- 0L
    Z1.S <- matrix(1, nrow = n, ncol = 1L)
    colnames(Z1.S) <- "Intercept"
  } else {
    p.S <- ncol(Z.S)
    Z1.S <- cbind(Intercept = 1, Z.S)
  }

  p1.X <- ncol(Z1.X)
  p1.S <- ncol(Z1.S)

  if (!is.null(prop.sd)) {
    if (length(prop.sd) == 1L) {
      if (!is.numeric(prop.sd) || !is.finite(prop.sd) || prop.sd <= 0) {
        stop("If 'prop.sd' is a scalar, it must be a single positive number.",
             call. = FALSE)
      }
    } else {
      if (!is.matrix(prop.sd) || !is.numeric(prop.sd)) {
        stop("Argument 'prop.sd' must be either a positive scalar or a numeric square matrix.",
             call. = FALSE)
      }
      if (nrow(prop.sd) != ncol(prop.sd)) {
        stop("If 'prop.sd' is a matrix, it must be square.", call. = FALSE)
      }
      if (any(!is.finite(prop.sd))) {
        stop("All entries of 'prop.sd' must be finite.", call. = FALSE)
      }
      if (!isTRUE(all.equal(prop.sd, t(prop.sd), tolerance = 1e-12))) {
        stop("If 'prop.sd' is a matrix, it must be symmetric.", call. = FALSE)
      }
    }
  }

  list(
    d = d,
    L = L,
    R = R,
    Z.X = Z.X,
    Z.S = Z.S,
    Z1.X = Z1.X,
    Z1.S = Z1.S,
    p.X = p.X,
    p.S = p.S,
    p1.X = p1.X,
    p1.S = p1.S,
    w = w,
    prop.sd = prop.sd,
    omitted = omitted
  )
}
