# Delete c file and open for writing
c_filename = 'src/skew_grid.c'
file.create(c_filename)
fp <- file(c_filename, open='at')
writeLines('#include <stdio.h>\n#include "skew_grid.h"\n', fp)

min_points = 4
max_points = 20
line <- sprintf("int min_grid_points = %d;\nint max_grid_points = %d;\n", min_points, max_points)
writeLines(line, fp);

for (K in seq(min_points, max_points)) {

  # Create basic vectors
  k = 0:K
  u = k/K; u[K+1] = (K-0.5)/K
  v = 1-(1-u)^2
  x = qnorm(0.5 + 0.5*v)
  xu = qnorm(0.5 + 0.5*u)
  c = 0.5/dnorm(x)
  cu = 0.5/dnorm(xu)
  f_v = 0.5/sqrt(1-v)
  f_v_inv = 1.0/f_v
  Kf_v2_inv = 1.0/(f_v * f_v * K)
  f_v_prime = 0.5 * f_v / (1-v)
  v_names <- c('u', 'v', 'x', 'xu', 'c', 'cu', 'f_v', 'f_v_inv', 'Kf_v2_inv', 'f_v_prime')

  # Create matrices
  m_names <- c('x_plus', 'x_minus', 'xu_plus', 'xu_minus')

  # Print lines for vectors
  for (v_name in v_names) {
    line <- sprintf("static double %s%d[%d] = {%s};", v_name, K, K+1,
                   paste(get(v_name), collapse=', '))
    writeLines(line, fp)
  }

  # Print lines for x_plus and x_minus matrices
  for (m_name in m_names) {
    writeLines(sprintf("static double %s%d[%d] = {", m_name, K, (K+1)*6), fp)
    for (i in seq(0, 5)) {
      s <- if (m_name=='x_plus' || m_name=='xu_plus') 1 else -1
      xvar <- if (m_name=='xu_plus' || m_name=='xu_minus') xu else x
      writeLines(sprintf("  %s,", paste((s*xvar)^i/factorial(i), collapse=', ')), fp)
    }
    writeLines("};", fp)
  }
  line <- sprintf("Grid g%d = {%d, %f, %f, %s, %s};\n", K, K, 1/K, log(K),
                  paste(sprintf("%s%d", v_names, K), collapse=', '),
                  paste(sprintf("%s%d", m_names, K), collapse=', '))
  writeLines(line, fp)
}

pointer_list = sprintf("%s, %s",
                       paste(rep("NULL", min_points), collapse=', '),
                       paste(sprintf("&g%d", min_points:max_points), collapse=', '))
line <- sprintf("Grid *grids[%d] = {%s};", max_points+1, pointer_list)
writeLines(line, fp)
close(fp)

