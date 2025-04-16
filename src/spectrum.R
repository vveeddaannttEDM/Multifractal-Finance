plot_multifractal_spectrum <- function(q, Hq) {
  # q: vector of q-values
  # Hq: corresponding Hurst exponents H(q)

  # Step 1: Compute τ(q)
  tq <- q * Hq - 1

  # Step 2: Compute α and f(α) using numerical derivatives
  dq <- diff(q)
  dtq <- diff(tq)
  alpha <- dtq / dq            # α(q) = dτ(q)/dq
  f_alpha <- q[-length(q)] * alpha - tq[-length(tq)]

  # Step 3: Plot the multifractal spectrum
  plot(alpha, f_alpha, type = "l", lwd = 2, col = "darkred",
       xlab = expression(alpha),
       ylab = expression(f(alpha)),
       main = "Multifractal Spectrum",
       xlim = range(alpha), ylim = range(f_alpha))
  grid()
}
