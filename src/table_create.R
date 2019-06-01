# Delete c file and open for writing
c_filename = 'src/skew_grid.c'
file.create(c_filename)
fp <- file(c_filename, open='at')
writeLines('#include "skew_grid.h"\n', fp)

for (K in seq(4, 20)) {

  # Create basic vectors
  k = 0:K
  u = k/K; u[K+1] = (K-0.5)/K
  v = 1-(1-u)^2
  x = qnorm(0.5 + 0.5*v)
  f_v = 0.5/sqrt(1-v)
  f_v_prime = 0.5 * f_v / (1-v)
  v_names <- c('u', 'v', 'x', 'f_v', 'f_v_prime')

  # Create matrices
  z_plus <- 0
  z_minus <- 0
  m_names <- c('z_plus', 'z_minus')

  # Print lines for vectors
  for (v_name in c('u', 'v', 'x', 'f_v', 'f_v_prime')) {
    line <- sprintf("static double %s%d[%d] = {%s};", v_name, K, K+1,
                   paste(get(v_name), collapse=', '))
    writeLines(line, fp)
  }

  # Print lines for z_plus and z_minus matrices
  for (m_name in m_names) {
    writeLines(sprintf("static double %s%d[%d] = {", m_name, K, (K+1)*6), fp)
    for (i in seq(0, 5)) {
      s <- if (m_name=='z_plus') 1 else -1
      writeLines(sprintf("  %s,", paste((s*x)^i, collapse=', ')), fp)
    }
    writeLines("};", fp)
  }
  line <- sprintf("Grid g%d = {%d, %s, %s};\n", K, K,
                  paste(sprintf("%s%d", v_names, K), collapse=', '),
                  paste(sprintf("%s%d", m_names, K), collapse=', '))
  writeLines(line, fp)
}
close(fp)

