#' Compute efficiency information for a skew draw evaluation
#'
#' \code{eff} does ...
#' @param lnf, log normalized density evaluated on a grid
#' @param true_lnf, true log (but unnormalized) density evaluated on the same grid
#' @param z_max, maximum grid point
#' @param n_pos, number of positive grid points (total number is 2*n_pos + 1)
#' @return A list of outputs
#' @export
eff <- function(lnf, true_lnf, z_max, n_pos) {
  f_true <- exp(true_lnf - max(true_lnf))
  f_true <- f_true / sum(f_true)
  integral <- sum(exp(lnf))*(z_max/n_pos)
  f <- exp(lnf - max(lnf))
  f <- f/sum(f)
  w <- f_true/f
  E_w <- sum(w * f_true)
  var_w <- sum((w-E_w)^2 * f_true)
  sd_w <- sqrt(var_w)
  skew_w <- sum((w-E_w)^3 * f_true) / sd_w^3
  list(integral=integral, E_w=E_w, var_w=var_w, sd_w=sd_w, skew_w=skew_w)
}
