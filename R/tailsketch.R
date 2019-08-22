ae_bar0 = -6.866007; ae_bar1 = -7.960650; ae_bar2 = -1.713753; ae_bar3 = 7.510672
ao_bar0 = -0.867413; ao_bar1 = -1.150060; ao_bar2 = 0.529299; ao_bar3 = 6.033791

x = seq(0, 2, by=0.01)
plot(x, ae_bar0 + ae_bar1*x + ae_bar2*x^2/2 + ae_bar3*x^3/6, ylim=c(-30, 0), type='l')
lines(x, ao_bar0 + ao_bar1*x + ao_bar2*x^2/2 + ao_bar3*x^3/6)
