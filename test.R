# Test discrete draw
d <- alias_cpp(c(1,2,3,4,5,4,3,2,1), 1000)
hist(d, breaks = seq(-0.5, 8.5, by=1))

# Test spline evaluation
u <- seq(0, 1, by=0.01)
f_u <- spline_eval_cpp(c(1, 2, 1), c(1, 0, -1), u)
plot(u, f_u, type='l')

# Test spline draw
u_draw <- spline_draw_cpp(c(1, 2, 1), c(1, 0, -1), 10000)
hist(u_draw)

# Test skew_eval_cpp
n <- 10
y_bar <- 5
theta <- 5
omega <- 16
x_max <- 1
r <- Inf; mu <- n*(theta-y_bar)/omega
#r <- 10; mu <- n*r*(theta-y_bar)/(omega*(r+theta))
case <- Po_GaPo(n, y_bar, r, theta, omega, mode=0, x_max, n_pos=1000)
plot(case$z, case$phi - 0.5*omega*case$z^2, type='l')
lines(case$z, case$phi, col='red')
lines(case$z, -0.5*case$x^2, col='green')

lnf <- skew_eval_cpp(0, case$h, mu, omega, case$z)
lines(case$z, lnf + 1.2, col='blue')
x_knots = c(0, 0.392830813649729, 0.764709673786387, 1.15034938037601, 1.59321881802305, 2.20041058121003, 2.69949670022498)
z_knots = x_knots / sqrt(omega)
abline(v=z_knots)
abline(v=-z_knots)

plot(case$z, lnf - case$phi + 0.5*case$x^2, type='l')
abline(v=z_knots)
abline(v=-z_knots)

