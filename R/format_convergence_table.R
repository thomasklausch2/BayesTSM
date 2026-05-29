#' Convergence table
#' @noRd
format_convergence_table <- function(g_x, g_s, eff_x, eff_s, digits_r = 3, digits_eff = 1) {
  tab_x <- data.frame(
    block = "X",
    parameter = names(g_x),
    R_hat = unname(round(g_x, digits_r)),
    ESS = unname(round(eff_x, digits_eff)),
    stringsAsFactors = FALSE
  )

  tab_s <- data.frame(
    block = "S",
    parameter = names(g_s),
    R_hat = unname(round(g_s, digits_r)),
    ESS = unname(round(eff_s, digits_eff)),
    stringsAsFactors = FALSE
  )

  out <- rbind(tab_x, tab_s)
  rownames(out) <- NULL
  out
}
