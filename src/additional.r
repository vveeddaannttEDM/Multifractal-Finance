# multifractal_analysis_extended.R
# ====================================================
# This script performs multifractal analysis of the financial volatility 
# (absolute returns) of the S&P 500 using a partition function approach.
#
# It extends the original analysis by adding:
#   1. Rolling window analysis to assess the time-evolution of scaling exponents.
#   2. Surrogate analysis (data shuffling) to test the significance of multifractality.
#   3. Bootstrap resampling to obtain confidence intervals for the scaling exponents.
#
# Before running, ensure you have the required packages installed:
# install.packages(c('tidyverse', 'pracma', 'boot'))
# Also, ensure that your data file (sp500_minutely.csv) is placed in the 'data' folder.

library(tidyverse)
library(pracma)  # For multifractal analysis functions
library(boot)

# ------------------------------
# 1. Load and Preprocess Data
# ------------------------------
# Load S&P 500 Minutely Data
# (Assume data has columns: Date, Close)
data <- read.csv('data/sp500_minutely.csv')
data$Date <- as.POSIXct(data$Date, format='%Y-%m-%d %H:%M:%S')

# Compute minutely returns and remove NA's
data <- data %>% 
  mutate(Return = log(Close / lag(Close))) %>% 
  na.omit()

# Define volatility as absolute returns
volatility <- abs(data$Return)
T_total <- length(volatility)

# ------------------------------
# 2. Direct Measure Analysis
# ------------------------------
# Define scales and moment orders q
scales <- floor(exp(seq(log(10), log(10000), length.out = 50)))
q_values <- seq(-3, 6, by = 1)

# Function to compute partition function for a given q and scale s:
partition_function <- function(q, s, vol) {
  N <- floor(length(vol) / s)
  # Compute the measure in each box:
  mu <- sapply(1:N, function(n) {
    sum(vol[((n - 1) * s + 1):(n * s)]) / sum(vol)
  })
  sum(mu^q)
}

# Compute the partition function for each q and each scale:
chi_q <- matrix(NA, nrow = length(q_values), ncol = length(scales))
for (i in seq_along(q_values)) {
  for (j in seq_along(scales)) {
    chi_q[i, j] <- partition_function(q_values[i], scales[j], volatility)
  }
}

# Estimate scaling exponents tau(q) by linear regression on log-log plots:
tau_q <- rep(NA, length(q_values))
for (i in seq_along(q_values)) {
  fit <- lm(log(chi_q[i, ]) ~ log(scales))
  tau_q[i] <- coef(fit)[2]
}

# Plot tau(q) vs q:
plot(q_values, tau_q, type = 'b', pch = 19, col = 'blue',
     main = 'Scaling Exponents tau(q)', xlab = 'q', ylab = 'tau(q)')

# ------------------------------
# 3. Inverse Measure Analysis (Exit Times)
# ------------------------------
# Define thresholds for exit time analysis:
thresholds <- seq(0.0001, 0.01, length.out = 50)

# Function to compute exit times given a threshold:
compute_exit_times <- function(vol, threshold) {
  exit_times <- c()
  cum_vol <- 0
  t <- 1
  while (t <= length(vol)) {
    cum_vol <- cum_vol + vol[t]
    if (cum_vol >= threshold) {
      exit_times <- c(exit_times, t)
      cum_vol <- 0
    }
    t <- t + 1
  }
  return(exit_times)
}

# Compute the inverse partition function for moment order p:
inverse_partition_function <- function(p, threshold, vol) {
  exit_times <- compute_exit_times(vol, threshold)
  if (length(exit_times) == 0) return(NA)
  mu_star <- exit_times / sum(exit_times)
  sum(mu_star^p)
}

p_values <- seq(-3, 6, by = 1)
chi_p_star <- matrix(NA, nrow = length(p_values), ncol = length(thresholds))
for (i in seq_along(p_values)) {
  for (j in seq_along(thresholds)) {
    chi_p_star[i, j] <- inverse_partition_function(p_values[i], thresholds[j], volatility)
  }
}

# Estimate scaling exponents theta(p) by linear regression:
theta_p <- rep(NA, length(p_values))
for (i in seq_along(p_values)) {
  valid <- !is.na(chi_p_star[i, ])
  fit <- lm(log(chi_p_star[i, valid]) ~ log(thresholds[valid]))
  theta_p[i] <- coef(fit)[2]
}

# Plot theta(p) vs p:
plot(p_values, theta_p, type = 'b', pch = 19, col = 'red',
     main = 'Scaling Exponents theta(p)', xlab = 'p', ylab = 'theta(p)')

# ------------------------------
# 4. Verify Inversion Formula
# ------------------------------
# The inversion formula predicts that tau(q) = -theta^{-1}(-q).
# We use interpolation to approximate theta^{-1}(-q):
theta_inv <- -approx(p_values, theta_p, xout = -q_values)$y

# Plot comparison:
plot(q_values, tau_q, type = 'b', col = 'blue', pch = 19,
     ylim = range(c(tau_q, theta_inv)),
     ylab = 'Exponents', main = 'Verification of Inversion Formula')
lines(q_values, theta_inv, type = 'b', col = 'green', pch = 17)
legend('topright', legend = c('tau(q)', '-theta^{-1}(-q)'),
       col = c('blue', 'green'), pch = c(19, 17))

# ------------------------------
# 5. Additional Suggestion 1: Rolling Window Analysis
# ------------------------------
# Break the volatility series into overlapping windows to see if tau(q) is stable over time.
window_size <- 10000  # e.g., 10,000 data points per window
step_size <- 5000     # overlapping windows (50% overlap)
num_windows <- floor((T_total - window_size) / step_size) + 1
tau_rolling <- matrix(NA, nrow = length(q_values), ncol = num_windows)
window_indices <- vector('list', num_windows)

for (w in 1:num_windows) {
  start_index <- ((w - 1) * step_size) + 1
  end_index <- start_index + window_size - 1
  window_indices[[w]] <- start_index:end_index
  vol_window <- volatility[start_index:end_index]
  
  for (i in seq_along(q_values)) {
    chi_temp <- sapply(scales, function(s) {
      if (floor(length(vol_window) / s) > 0) {
        partition_function(q_values[i], s, vol_window)
      } else {
        NA
      }
    })
    valid <- !is.na(chi_temp)
    if (sum(valid) > 2) {
      fit <- lm(log(chi_temp[valid]) ~ log(scales[valid]))
      tau_rolling[i, w] <- coef(fit)[2]
    } else {
      tau_rolling[i, w] <- NA
    }
  }
}

# For example, plot the evolution of tau(2) over time:
q_index <- which(q_values == 2)
time_points <- sapply(window_indices, function(idx) mean(idx))
plot(time_points, tau_rolling[q_index, ], type = 'b', pch = 19, col = 'purple',
     main = 'Rolling Window Analysis: tau(2) over Time',
     xlab = 'Time Index', ylab = 'tau(2)')

# ------------------------------
# 6. Additional Suggestion 2: Surrogate Analysis
# ------------------------------
# Shuffle the volatility series to remove temporal correlations and recalc tau(q).
set.seed(123)  # For reproducibility
volatility_shuffled <- sample(volatility)
chi_q_surrogate <- matrix(NA, nrow = length(q_values), ncol = length(scales))
for (i in seq_along(q_values)) {
  for (j in seq_along(scales)) {
    chi_q_surrogate[i, j] <- partition_function(q_values[i], scales[j], volatility_shuffled)
  }
}
tau_q_surrogate <- rep(NA, length(q_values))
for (i in seq_along(q_values)) {
  fit <- lm(log(chi_q_surrogate[i, ]) ~ log(scales))
  tau_q_surrogate[i] <- coef(fit)[2]
}

# Plot comparison: Original vs. Surrogate tau(q)
plot(q_values, tau_q, type = 'b', pch = 19, col = 'blue',
     main = 'Original vs. Surrogate tau(q)', xlab = 'q', ylab = 'tau(q)')
lines(q_values, tau_q_surrogate, type = 'b', pch = 17, col = 'orange')
legend('topright', legend = c('Original', 'Surrogate'), col = c('blue', 'orange'), pch = c(19, 17))

# ------------------------------
# 7. Additional Suggestion 3: Bootstrap Resampling for Confidence Intervals
# ------------------------------
# Define a function for bootstrapping tau(q) for a given moment q.
bootstrap_tau <- function(data, indices, q, scales) {
  vol_sample <- data[indices]
  chi_vals <- sapply(scales, function(s) partition_function(q, s, vol_sample))
  fit <- lm(log(chi_vals) ~ log(scales))
  coef(fit)[2]
}

B <- 100  # Number of bootstrap replications
tau_bootstrap <- matrix(NA, nrow = length(q_values), ncol = B)
for (i in seq_along(q_values)) {
  boot_result <- boot(data = volatility, statistic = bootstrap_tau,
                      R = B, q = q_values[i], scales = scales)
  tau_bootstrap[i, ] <- boot_result$t
}

# Compute 95% confidence intervals for each q:
tau_ci <- t(apply(tau_bootstrap, 1, quantile, probs = c(0.025, 0.975)))

# Plot tau(q) with error bars:
plot(q_values, tau_q, type = 'b', pch = 19, col = 'blue',
     ylim = range(c(tau_q, tau_ci)),
     main = 'tau(q) with Bootstrap 95% Confidence Intervals',
     xlab = 'q', ylab = 'tau(q)')
for (i in seq_along(q_values)) {
  arrows(q_values[i], tau_ci[i, 1], q_values[i], tau_ci[i, 2],
         code = 3, angle = 90, length = 0.05, col = 'red')
}

# ------------------------------
# 8. Save Results (Optional)
# ------------------------------
# Save plots to a PNG file in the 'results' folder:
png('results/scaling_exponents_plots.png', width = 800, height = 800)
par(mfrow = c(2, 2))
plot(q_values, tau_q, type = 'b', pch = 19, col = 'blue',
     main = 'tau(q)', xlab = 'q', ylab = 'tau(q)')
plot(p_values, theta_p, type = 'b', pch = 19, col = 'red',
     main = 'theta(p)', xlab = 'p', ylab = 'theta(p)')
plot(q_values, tau_q, type = 'b', col = 'blue', pch = 19,
     ylim = range(c(tau_q, theta_inv)), ylab = 'Exponents',
     main = 'tau(q) vs. -theta^{-1}(-q)')
lines(q_values, theta_inv, type = 'b', col = 'green', pch = 17)
plot(time_points, tau_rolling[q_index, ], type = 'b', pch = 19, col = 'purple',
     main = 'Rolling Window: tau(2)', xlab = 'Time Index', ylab = 'tau(2)')
dev.off()
