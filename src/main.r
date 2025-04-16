# Required Libraries
library(tidyverse)
library(pracma)  # For multifractal analysis functions

# Load S&P 500 Minutely Data (replace 'sp500_minutely.csv' with your dataset)
data <- read.csv('sp500_minutely.csv')
data$Date <- as.POSIXct(data$Date, format='%Y-%m-%d %H:%M:%S')

# Compute Minutely Returns
data <- data %>% mutate(Return = log(Close / lag(Close))) %>% na.omit()

# Direct Volatility Measure
volatility <- abs(data$Return)
T <- length(volatility)

# Define scales for partition function
scales <- floor(exp(seq(log(10), log(10000), length.out = 50)))

# Compute Partition Function for Direct Measure
partition_function <- function(q, s, volatility) {
  N <- floor(T / s)
  mu <- sapply(1:N, function(n) sum(volatility[((n-1)*s + 1):(n*s)]) / sum(volatility))
  sum(mu^q)
}

# Calculate χ_q(s) for multiple q values
q_values <- seq(-3, 6, by=1)
chi_q <- sapply(q_values, function(q) sapply(scales, function(s) partition_function(q, s, volatility)))

# Estimate Scaling Exponents τ(q)
tau_q <- apply(chi_q, 1, function(chi) {
  lm(log(chi) ~ log(scales))$coefficients[2]  # Slope gives τ(q)
})

# Plot τ(q) vs q
plot(q_values, tau_q, type='b', pch=19, col='blue', main='Scaling Exponents τ(q)', xlab='q', ylab='τ(q)')

# Inverse Measure (Exit Time Analysis)
thresholds <- seq(0.0001, 0.01, length.out = 50)

# Compute Exit Times
compute_exit_times <- function(volatility, threshold) {
  exit_times <- c()
  cum_vol <- 0
  t <- 1
  while (t <= length(volatility)) {
    cum_vol <- cum_vol + volatility[t]
    if (cum_vol >= threshold) {
      exit_times <- c(exit_times, t)
      cum_vol <- 0
    }
    t <- t + 1
  }
  exit_times
}

# Compute Partition Function for Inverse Measure
inverse_partition_function <- function(p, threshold, volatility) {
  exit_times <- compute_exit_times(volatility, threshold)
  mu_star <- exit_times / sum(exit_times)
  sum(mu_star^p)
}

# Calculate χ*_p(Δv) for multiple p values
p_values <- seq(-3, 6, by=1)
chi_p_star <- sapply(p_values, function(p) sapply(thresholds, function(th) inverse_partition_function(p, th, volatility)))

# Estimate Scaling Exponents θ(p)
theta_p <- apply(chi_p_star, 1, function(chi_star) {
  lm(log(chi_star) ~ log(thresholds))$coefficients[2]  # Slope gives θ(p)
})

# Plot θ(p) vs p
plot(p_values, theta_p, type='b', pch=19, col='red', main='Scaling Exponents θ(p)', xlab='p', ylab='θ(p)')

# Verify Inversion Formula
tau_inversion <- -approx(p_values, theta_p, xout = -q_values)$y

# Plotting τ(q) and -θ^-1(-q)
plot(q_values, tau_q, type='b', col='blue', pch=19, ylim=range(c(tau_q, tau_inversion)), ylab='Exponents', main='Inversion Formula Verification')
lines(q_values, tau_inversion, type='b', col='green', pch=17)
legend('topright', legend=c('τ(q)', '-θ^{-1}(-q)'), col=c('blue', 'green'), pch=c(19,17))
