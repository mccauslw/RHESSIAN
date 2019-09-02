library(tidyverse)
set.seed(5750418)
n_cases = 500
tbl <- tibble(r = rgamma(n_cases, 2, rate=0.1),
              n = 1+rnbinom(n_cases, size=10, mu=3),
              theta = rgamma(n_cases, 2, rate=1),
              y_bar = rnbinom(n_cases, size=r/4, mu=theta)/n,
              omega = rgamma(n_cases, 2, rate=0.2),
              a2 = 0, a3 = 0, a4 = 0, a5 = 0,
              e_old = 0, e20_0 = 0, e20_1 = 0)

#Po_GaPo <- function(n, y_bar, r, theta, omega, mode=0, x_max, n_pos=1000)

#  list(h = c(0, 0, h_2, h_3, h_4, h_5), a = a, x = x, z = z,
#       phi = phi, phi_m = phi_m, phi_o = phi_o, phi_e = phi_e,
#       phi_Taylor = phi_Taylor, phi_e_Taylor = phi_e_Taylor)

x_max <- 5.0
n_pos <- 500
mode <- 0
for (i in 1:n_cases) {

  n = tbl$n[i]; y_bar = tbl$y_bar[i]; theta = tbl$theta[i]; omega = tbl$omega[i]; r = tbl$r[i]
  mu <- n*r*(theta-y_bar)/(omega*(r+theta))
  z_max <- x_max / sqrt(omega)
  case_info <- Po_GaPo(n, y_bar, r, theta, omega, mode, x_max, n_pos)
  tbl$a2[i] <- case_info$a[3]
  tbl$a3[i] <- case_info$a[4]
  tbl$a4[i] <- case_info$a[5]
  tbl$a5[i] <- case_info$a[6]
  true_lnf = case_info$phi - 0.5*omega*case_info$z^2

  lnf <- skew_eval_cpp(K, 2, mode, case_info$h, mu, omega, case_info$z)
  eff_results <- eff(lnf, true_lnf, z_max, n_pos)
  tbl$e_old[i] <- eff_results$sd_w
  for (K in c(20)) {
    for (code in c(0, 1)) {
      name = sprintf("e%d_%d", K, code)
      lnf <- skew_eval_cpp(K, code, mode, case_info$h, mu, omega, case_info$z)
      eff_results <- eff(lnf, true_lnf, z_max, n_pos)
      tbl[[name]][i] <- eff_results$sd_w
    }
  }
}

tbl_screen = filter(tbl, omega < 5.0)
quantile(tbl_screen$e20_1, seq(0.1, 0.9, by=0.1))

tbl <- mutate(tbl, ln_adv = log(e_old) - log(e20_1),
              a3n = a3/(-a2)^1.5, a4n = a4/(-a2)^2, a5n = a5/(-a2)^2.5,
              e20_0_best = e20_0 < pmin(e20_1, e_old),
              e20_1_best = e20_1 < pmin(e20_0, e_old))

mylm <- lm(formula = ln_adv ~ log(-a2) + abs(a3n) + abs(a4n), data = tbl)
summary(mylm)

plot(tbl$a2, tbl$a3, pch=20, xlim=c(-3, 0), ylim = c(-1, 0.5))
points(tbl$a2[tbl$e20_1_best], tbl$a3[tbl$e20_1_best], col='red')
points(tbl$a2[tbl$e20_0_best], tbl$a3[tbl$e20_0_best], col='green')


