# Test skew_eval_cpp
n <- 1
y_bar <- 1
theta <- 1
omega <- 100
x_max <- 5.0; z_max <- x_max / sqrt(omega)
n_pos <- 20000
#r <- 20; mu <- n*(theta-y_bar)/omega
r <- 10; mu <- n*r*(theta-y_bar)/(omega*(r+theta))
case <- Po_GaPo(n, y_bar, r, theta, omega, mode=0, x_max, n_pos=n_pos)
true_lnf = case$phi - 0.5*omega*case$z^2
plot(case$z, true_lnf, type='l')
lines(case$z, case$phi, col='red')
lines(case$z, -0.5*case$x^2, col='green')
lines(case$z, case$phi_Taylor, col='purple')

K <- 8
code = 1
mode = 0

print(system.time(for(i in 1:100) li_new <- skew_draw_cpp(K, code, mode, case$h, mu, omega, 40000)))
print(system.time(for(i in 1:100) li_old <- skew_draw_cpp(K, 2, mode, case$h, mu, omega, 40000)))

print(system.time(for(i in 1:100) lnf_new <- skew_eval_cpp(K, code, mode, case$h, mu, omega, case$z)))
print(system.time(for(i in 1:100) lnf_old <- skew_eval_cpp(K, 2, mode, case$h, mu, omega, case$z)))

k <- 0:K
u <- k/K; u[K+1] = (K-0.5)/K
v <- 1-(1-u)^2
x_knots <- if (code==1) qnorm(0.5 + 0.5*v) else qnorm(0.5 + 0.5*u)
z_knots <- x_knots / sqrt(omega)
lines(case$z, lnf_new, col='red')
abline(v=z_knots)
abline(v=-z_knots)

plot(case$z, lnf_new - true_lnf, type='l', col='red')
lines(case$z, lnf_old - true_lnf, type='l', col='green')
abline(v=z_knots)
abline(v=-z_knots)

hist(li_new$draws, 100, freq=FALSE)
lines(case$z, exp(lnf_new))
abline(v=z_knots)
abline(v=-z_knots)

for (i in 1:2) {
  lnf <- if (i==1) lnf_new else lnf_old
  print(sum(exp(lnf))*(z_max/n_pos))
  diff = true_lnf[n_pos+1] - lnf[n_pos+1]
  Ew = sum(exp(true_lnf-lnf-diff) * exp(lnf)*(z_max/n_pos))
  Ew2 = sum(exp(2*(true_lnf-lnf-diff)) * exp(lnf)*(z_max/n_pos))
  var_w = Ew2 - Ew^2
  sd_w = sqrt(var_w)
  print(sprintf("Ew: %f, Ew2: %f, sd_w: %f", Ew, Ew2, sd_w))
}
