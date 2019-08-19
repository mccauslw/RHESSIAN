# Test skew_eval_cpp
n <- 2
y_bar <- 0.5
theta <- 2.92
omega <- 1.09
#r <- 20; mu <- n*(theta-y_bar)/omega
r <- 16.4; mu <- n*r*(theta-y_bar)/(omega*(r+theta))

x_max <- 5.0; z_max <- x_max / sqrt(omega)
n_pos <- 200 # 20000
case <- Po_GaPo(n, y_bar, r, theta, omega, mode=0, x_max, n_pos=n_pos)
true_lnf = case$phi - 0.5*omega*case$z^2
plot(case$z, true_lnf, type='l', ylim=c(-30,0))
#lines(case$z, case$phi, col='red')
#lines(case$z, -0.5*case$x^2, col='green')
lines(case$z, case$phi_Taylor, col='purple')

K <- 16
code = 1
mode = 0

print(system.time(for(i in 1:100) li_new <- skew_draw_cpp(K, code, mode, case$h, mu, omega, 40000)))
print(system.time(for(i in 1:100) li_old <- skew_draw_cpp(K, 2, mode, case$h, mu, omega, 40000)))

print(system.time(for(i in 1:100) new_lnf <- skew_eval_cpp(K, code, mode, case$h, mu, omega, case$z)))
print(system.time(for(i in 1:100) old_lnf <- skew_eval_cpp(K, 2, mode, case$h, mu, omega, case$z)))

k <- 0:K
u <- k/K; u[K+1] = (K-0.5)/K
v <- 1-(1-u)^2
x_knots <- if (code==1) qnorm(0.5 + 0.5*v) else qnorm(0.5 + 0.5*u)
z_knots <- x_knots / sqrt(omega)
lines(case$z, new_lnf, col='red')
abline(v=z_knots)
abline(v=-z_knots)

plot(case$z, new_lnf - true_lnf, type='l', col='red', xlim=c(-z_max, z_max), ylim=c(-2,2))
lines(case$z, old_lnf - true_lnf, type='l', col='green')
abline(v=z_knots)
abline(v=-z_knots)

hist(li_new$draws, 100, freq=FALSE)
lines(case$z, exp(new_lnf))
points(li_new$draws[1:100], exp(li_new$ln_f[1:100]), pch=20, col='red')
abline(v=z_knots)
abline(v=-z_knots)

print(eff(new_lnf, true_lnf, z_max, n_pos))
print(eff(old_lnf, true_lnf, z_max, n_pos))
