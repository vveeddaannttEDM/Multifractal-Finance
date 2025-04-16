hurst_rs <- function(series, min_window = 10, max_window = 100, step = 10) {
  N <- length(series)
  hurst_data <- data.frame()

  for (window in seq(min_window, max_window, by = step)) {
    n_segments <- floor(N / window)
    rs_values <- numeric(n_segments)

    for (i in 1:n_segments) {
      segment <- series[((i - 1) * window + 1):(i * window)]
      segment <- segment - mean(segment)
      Z <- cumsum(segment)
      R <- max(Z) - min(Z)
      S <- sd(segment)

      if (S != 0) {
        rs_values[i] <- R / S
      } else {
        rs_values[i] <- NA
      }
    }

    rs_mean <- mean(rs_values, na.rm = TRUE)
    hurst_data <- rbind(hurst_data, data.frame(window = window, RS = rs_mean))
  }

  # Linear regression on log-log plot
  log_rs <- log(hurst_data$RS)
  log_window <- log(hurst_data$window)

  model <- lm(log_rs ~ log_window)
  hurst_exp <- coef(model)[2]

  plot(log_window, log_rs, main = "R/S Plot", xlab = "log(Window)", ylab = "log(R/S)", pch = 20)
  abline(model, col = "blue", lwd = 2)

  return(hurst_exp)
}
