library(tidyverse)
set.seed(5750418)
n_cases = 200
tbl <- tibble(r = rgamma(n_cases, 2, rate=0.1),
              n = 1+rnbinom(n_cases, size=10, mu=3),
              theta = rgamma(n_cases, 2, rate=1),
              y_bar = rnbinom(n_cases, size=r/4, mu=theta)/n,
              omega = rgamma(n_cases, 2, rate=0.2),
              a2 = 0, a3 = 0, a4 = 0, a5 = 0,
              e_old = 0, e4_0 = 0, e4_1 = 0, e8_0 = 0, e8_1 = 0, e16_0 = 0, e16_1 = 0)

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
  #print(eff_results)
  tbl$e_old[i] <- eff_results$sd_w
  for (K in c(4, 8, 16)) {
    for (code in c(0, 1)) {
      name = sprintf("e%d_%d", K, code)
      lnf <- skew_eval_cpp(K, code, mode, case_info$h, mu, omega, case_info$z)
      eff_results <- eff(lnf, true_lnf, z_max, n_pos)
      #print(eff_results)
      tbl[[name]][i] <- eff_results$sd_w
    }
  }
}

tbl_screen <- filter(tbl, e_old<1)

