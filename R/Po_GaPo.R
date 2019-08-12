#' Compute info for Poisson and Gamma-Poisson models
#'
#' \code{Po_GoPo} does ...
#' @param n, number of observations
#' @param y_bar, mean count
#' @param r, Gamma-Poisson parameter, r=inf is Poisson special case
#' @param theta, multiplicative scale parameter for Poisson or Gamma-Poisson
#' @param mode, mode of distribution
#' @param omega, prior precision
#' @param x_max, maximum grid-point value
#' @param n_pos, number of positive grid points
#' @return A list of outputs
#' @export
Po_GaPo <- function(n, y_bar, r, theta, omega, mode=0, x_max, n_pos=1000) {
  # The prior is x ~ N(mu, sigma^2). Let omega = 1/sigma^2
  # The sample is Y=(Y_1,...,Y_n), with, for i=1,...,n,
  #   Y_i ~ iid Po(theta e^x)
  #   or Y_i ~ iid GaPo(r, (theta/r) e^x),
  # Note that the Poisson case is the limit of the Gamma Poisson case as r->infty.
  #
  # In both cases, the sufficient statistic is (n, y_bar), where y_bar is the
  # sample mean.
  # The parameter values omega, theta and r, and the sufficient stastistic
  # (n, y_bar) are specified here below.
  # The value of mu is set to make the posterior mode of x equal to zero.
  sigma = omega^-0.5

  # Step two: Compute derivatives of phi = log f(x|y) - log f_N(0, omega^{-1}),
  # normalized so that phi(0) = 0
  if (is.infinite(r)) {
    # Value of mu parameter giving a mode at x=0
    mu = n*(theta-y_bar)/omega
    # Derivatives of phi at zero, from 2nd to 5th
    h_2 = -n*theta
    h_3 = -n*theta
    h_4 = -n*theta
    h_5 = -n*theta
  } else {
    # Value of theta parameter giving a mode at x=0
    mu = n*r*(theta-y_bar)/(omega*(r+theta))
    th_by_r = theta/r
    # Derivatives of phi at zero, from 2nd to 5th
    h_2 = -n*(r+y_bar)*th_by_r / (1+th_by_r)^2
    h_3 = -n*(r+y_bar)*(th_by_r - th_by_r^2) / (1+th_by_r)^3
    h_4 = -n*(r+y_bar)*(th_by_r - 4*th_by_r^2 + th_by_r^3) / (1+th_by_r)^4
    h_5 = -n*(r+y_bar)*(th_by_r - 11*th_by_r^2 + 11*th_by_r^3 - th_by_r^4) / (1+th_by_r)^5
  }

  # Normalized values of h_2
  a = c(0, 0, h_2/omega, h_3/omega^{1.5}, h_4/omega^2, h_5/omega^{2.5})

  # Step three: Compute true phi and its odd and even parts on a fine grid, as an
  # illustration.
  c = n*y_bar + omega*mu
  x = seq(-x_max, x_max, length.out = 2*n_pos+1)
  x.plus = x[x>=0]
  z = x*sigma
  if (is.infinite(r)) {
    phi = c*z - n*theta*(exp(z)-1)
    phi_m = -c*z - n*theta*(exp(-z)-1) # phi(-x)
    phi_o = c*z - n*theta*sinh(z)      # odd part of phi
    phi_e = n*theta*(1-cosh(z))        # even part of phi
  } else {
    phi = c*z - n*(r+y_bar)*log((1 + th_by_r * exp(z))/(1 + th_by_r))
    phi_m = -c*z - n*(r+y_bar)*log((1 + th_by_r * exp(-z))/(1 + th_by_r))
    phi_o = c*z - 0.5*n*(r+y_bar)*log((1 + th_by_r * exp(z))/(1 + th_by_r * exp(-z)))
    phi_e = -0.5*n*(r+y_bar)*(log(1 + th_by_r^2 + 2*th_by_r*cosh(z)) - 2*log(1 + th_by_r))
  }
  phi_e_Taylor = 0.5 * (h_2 - omega) * z^2 - h_4 * z^4/24 + log(cosh(h_3 * z^3/6 + h_5 * z^5/120))
  phi_Taylor = (h_2 - omega) * z^2/2 + h_3 * z^3/6 + h_4 * z^4/24 + h_5 * z^5/120

  list(h = c(0, 0, h_2, h_3, h_4, h_5), a = a, x = x, z = z,
       phi = phi, phi_m = phi_m, phi_o = phi_o, phi_e = phi_e,
       phi_Taylor = phi_Taylor, phi_e_Taylor = phi_e_Taylor)
}
