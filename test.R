# Test discrete draw
d = alias_cpp(c(1,2,3,4,5,4,3,2,1), 1000)
hist(d, breaks = seq(-0.5, 8.5, by=1))

# Test spline evaluation
u = seq(0, 1, by=0.01)
f_u = spline_eval_cpp(c(1, 2, 1), c(1, 0, -1), u)
plot(u, f_u, type='l')

# Test spline draw
u_draw = spline_draw_cpp(c(1, 2, 1), c(1, 0, -1), 10000)
hist(u_draw)
