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

t <- seq(0, 1, by=0.01)
h00 <- (1+2*t)*(1-t)^2
h10 <- t*(1-t)^2
h01 <- t^2*(3-2*t)
h11 <- t^2*(t-1)
plot(t, 345.6*h00 -1900.8*h10)

# Test skew_eval_cpp
n <- 5
y_bar <- 2
theta <- 2
omega <- 9
x_max <- 4.0; z_max <- x_max / sqrt(omega)
n_pos <- 200
#r <- 20; mu <- n*(theta-y_bar)/omega
r <- 10; mu <- n*r*(theta-y_bar)/(omega*(r+theta))
case <- Po_GaPo(n, y_bar, r, theta, omega, mode=0, x_max, n_pos=n_pos)
plot(case$z, case$phi - 0.5*omega*case$z^2, type='l')
lines(case$z, case$phi, col='red')
lines(case$z, -0.5*case$x^2, col='green')
lines(case$z, case$phi_Taylor, col='purple')

K <- 10
code = 2
mode = 0
lnf <- skew_eval_cpp(K, code, mode, case$h, mu, omega, case$z)

k <- 0:K
u <- k/K; u[K+1] = (K-0.5)/K
v <- 1-(1-u)^2
x_knots <- if (is_v) qnorm(0.5 + 0.5*v) else qnorm(0.5 + 0.5*u)
z_knots <- x_knots / sqrt(omega)
lines(case$z, lnf, col='blue')
abline(v=z_knots)
abline(v=-z_knots)

plot(case$z, lnf - case$phi + 0.5*case$x^2, type='l')
abline(v=z_knots)
abline(v=-z_knots)

print(sum(exp(lnf))*(z_max/n_pos))
