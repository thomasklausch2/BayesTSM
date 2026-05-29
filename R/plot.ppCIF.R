#' Plot posterior predictive CIFs
#'
#' Plot method for objects returned by \code{ppCIF()}.
#'
#' Produces posterior predictive CIF plots for \code{X}, \code{S}, and
#' \code{Y = X + S}. Regardless of whether the object was created with
#' \code{type = "quantiles"} or \code{type = "percentiles"}, the plot is shown
#' in CDF form with time on the x-axis and probability on the y-axis.
#'
#' If \code{x$type = "quantiles"}, the stored posterior predictive CDF values
#' are plotted directly. If \code{x$type = "percentiles"}, the stored posterior
#' predictive quantiles are plotted against the probability grid, so that the
#' resulting graph is again displayed in CDF form.
#'
#' @param x An object returned by \code{ppCIF()}.
#' @param y Ignored.
#' @param ci Logical. If \code{TRUE}, pointwise 95 percent credible regions
#'   are added as shaded areas.
#' @param layout Character string. Either \code{"separate"} for three
#'   horizontally adjacent panels, or \code{"combined"} for a single chart with
#'   all three curves.
#' @param main Optional plot title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param xlim Optional numeric vector of length two giving the x-axis limits.
#'   If \code{NULL}, the limits are chosen automatically.
#' @param facet_scales Character string passed to
#'   \code{ggplot2::facet_wrap(scales = ...)} when
#'   \code{layout = "separate"}. The default is \code{"free_x"}.
#' @param ... Further graphical arguments, currently ignored.
#'
#' @return Invisibly returns the resulting \code{ggplot} object.
#'
#' @seealso [bayestsm()] [ppCIF()]
#'
#' @examples
#' \dontrun{
#' # Obtain quantiles for provided percentiles
#' # For a bayestsm model output saved in mod_slice (see ?bayestsm examples), run:
#' postCDF_perc = ppCIF(mod_slice, warmup = 500, perc = seq(0, 0.99, 0.01),
#'                      method = 'simulation', type = 'percentiles')
#'
#' # straight forward plot over provided percentiles - quantiles combintation
#' plot(postCDF_perc, xlim = c(0,40))
#' }
#'
#' @method plot ppCIF
#' @export
plot.ppCIF <- function(x, y = NULL, ci = TRUE,
                       layout = c("separate", "combined"),
                       main = NULL, xlab = NULL, ylab = NULL,
                       xlim = NULL,
                       facet_scales = "free_x", ...) {

  layout <- match.arg(layout)

  if (!inherits(x, "ppCIF")) {
    stop("x must be an object of class 'ppCIF'.")
  }
  if (is.null(x$type)) {
    stop("The ppCIF object does not contain a 'type' element.")
  }
  if (is.null(x$grid)) {
    stop("The ppCIF object does not contain a 'grid' element.")
  }
  if (!is.null(xlim) && (!is.numeric(xlim) || length(xlim) != 2L || anyNA(xlim))) {
    stop("xlim must be NULL or a numeric vector of length two without missing values.")
  }

  make_polygon_df <- function(time_lower, time_upper, prob, outcome) {
    data.frame(
      time = c(time_lower, rev(time_upper)),
      prob = c(prob, rev(prob)),
      outcome = outcome
    )
  }

  if (x$type == "quantiles") {
    df <- rbind(
      data.frame(
        time = x$grid,
        prob = x$med.p.x,
        lower = x$p.x.ci[1, ],
        upper = x$p.x.ci[2, ],
        outcome = "X"
      ),
      data.frame(
        time = x$grid,
        prob = x$med.p.s,
        lower = x$p.s.ci[1, ],
        upper = x$p.s.ci[2, ],
        outcome = "S"
      ),
      data.frame(
        time = x$grid,
        prob = x$med.p.y,
        lower = x$p.y.ci[1, ],
        upper = x$p.y.ci[2, ],
        outcome = "Y = X + S"
      )
    )
    poly_df <- NULL

  } else if (x$type == "percentiles") {
    df <- rbind(
      data.frame(
        time = x$med.q.x,
        prob = x$grid,
        lower = x$q.x.ci[1, ],
        upper = x$q.x.ci[2, ],
        outcome = "X"
      ),
      data.frame(
        time = x$med.q.s,
        prob = x$grid,
        lower = x$q.s.ci[1, ],
        upper = x$q.s.ci[2, ],
        outcome = "S"
      ),
      data.frame(
        time = x$med.q.y,
        prob = x$grid,
        lower = x$q.y.ci[1, ],
        upper = x$q.y.ci[2, ],
        outcome = "Y = X + S"
      )
    )

    poly_df <- rbind(
      make_polygon_df(
        time_lower = x$q.x.ci[1, ],
        time_upper = x$q.x.ci[2, ],
        prob = x$grid,
        outcome = "X"
      ),
      make_polygon_df(
        time_lower = x$q.s.ci[1, ],
        time_upper = x$q.s.ci[2, ],
        prob = x$grid,
        outcome = "S"
      ),
      make_polygon_df(
        time_lower = x$q.y.ci[1, ],
        time_upper = x$q.y.ci[2, ],
        prob = x$grid,
        outcome = "Y = X + S"
      )
    )
  } else {
    stop("Unknown ppCIF type: ", x$type)
  }

  df$outcome <- factor(df$outcome, levels = c("X", "S", "Y = X + S"))
  if (!is.null(poly_df)) {
    poly_df$outcome <- factor(poly_df$outcome, levels = c("X", "S", "Y = X + S"))
  }

  if (is.null(xlab)) xlab <- "Time"
  if (is.null(ylab)) ylab <- "Probability"
  if (is.null(main)) main <- "Posterior predictive CIFs"

  if (layout == "separate") {

    p <- ggplot2::ggplot()

    if (ci) {
      if (x$type == "quantiles") {
        p <- p +
          ggplot2::geom_ribbon(
            data = df,
            ggplot2::aes(x = time, ymin = .data$lower, ymax = .data$upper),
            alpha = 0.2
          )
      } else {
        p <- p +
          ggplot2::geom_polygon(
            data = poly_df,
            ggplot2::aes(x = time, y = .data$prob, group = .data$outcome),
            alpha = 0.2
          )
      }
    }

    p <- p +
      ggplot2::geom_path(
        data = df,
        ggplot2::aes(x = time, y = .data$prob),
        linewidth = 0.8
      ) +
      ggplot2::facet_wrap(~ outcome, nrow = 1, scales = facet_scales) +
      ggplot2::labs(
        title = main,
        x = xlab,
        y = ylab
      ) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold")
      )

  } else {

    p <- ggplot2::ggplot()

    if (ci) {
      if (x$type == "quantiles") {
        p <- p +
          ggplot2::geom_ribbon(
            data = df,
            ggplot2::aes(x = time, ymin = .data$lower, ymax = .data$upper, fill = .data$outcome),
            alpha = 0.2,
            inherit.aes = FALSE
          )
      } else {
        p <- p +
          ggplot2::geom_polygon(
            data = poly_df,
            ggplot2::aes(x = time, y = .data$prob, group = .data$outcome, fill = .data$outcome),
            alpha = 0.2,
            inherit.aes = FALSE
          )
      }
    }

    p <- p +
      ggplot2::geom_path(
        data = df,
        ggplot2::aes(
          x = time,
          y = .data$prob,
          colour = .data$outcome,
          linetype = .data$outcome,
          group = .data$outcome
        ),
        linewidth = 0.9
      ) +
      ggplot2::labs(
        title = main,
        x = xlab,
        y = ylab,
        colour = "Outcome",
        linetype = "Outcome",
        fill = "Outcome"
      ) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1)) +
      ggplot2::theme_bw()
  }

  print(p)
  invisible(p)
}
