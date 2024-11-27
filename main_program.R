library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(ggcorrplot)
library(MASS)
library(stringr)
library(car)
library(olsrr)
library(stabledist)
library(readxl)
library(stats)
library(forecast)
library(zoo)
library(tseries)
library(gridExtra)
library(LSMonteCarlo)

####################################### GBM ####################################
# set.seed(21)
S0 <- 100 # Underlying asset price
r <- 0.05 # risk-free rate
nsteps <- 12 # Time steps
npaths <- 2000 # Number of paths
T <- 1 # Maturity (years)
dt <- T/nsteps # Discritization (time)
K <- 90 # Strike price
sigma <- 0.2 # Volatility
discount <- exp(-r * dt) # For discounted payoff
#####

# Simulating stock price
dW <- matrix(rnorm(npaths * nsteps, mean = 0, sd = sqrt(dt)), ncol = nsteps)
S <- matrix(0, nrow = npaths, ncol = nsteps + 1)

S[, 1] <- S0
for (t in 1:nsteps) {
  S[, t + 1] <- S[, t] * exp((r - 0.5 * sigma^2) * dt + sigma * dW[, t])
}

# Payoff matrix for a Bermudan Call option
payoff <- pmax(S - K, 0)
# Initialize cashflow matrix
cashflow <- matrix(nrow = npaths, ncol = (nsteps + 1))

############################# Longstaff Approach ###############################
# Specify exercise frequency (Change "by" argument for time jumps)
exercise_points <- seq(1, nsteps, by = 1)

# Backward induction to compute continuation values (ORG)
# t follows col numbers NOT time. So for discounts t is changed to match the time
# of the column. (There are 13 columns and 12 times. Col 1 represents time 0).
for (t in (nsteps + 1):1) {
  # Skip steps that are not exercise points for Bermudan options
  if (!(t %in% exercise_points)) {
    
    cashflow[, t] <- payoff[, t] * discount^(t-1)
    next
  }
  if (t == (nsteps + 1)) {
    
    cashflow[, t] <- payoff[, t]
    next
  }
  # Only consider in-the-money paths
  in_the_money <- which(payoff[, t] > 0)
  if (length(in_the_money) > 0) {
    
  # Fit regression to estimate continuation value
  X <- S[in_the_money, t]
  Y <- payoff[in_the_money, t + 1] * discount^(t)
  continuation_value <- predict(lm(Y ~ I(X) + I(X^2)),
                                newdata = data.frame(X = S[in_the_money, t]))
  
  # Exercise or continue decision
  exercise <- payoff[in_the_money, t] > continuation_value
  cashflow[exercise, t] <- payoff[exercise, t]
  cashflow[exercise, (t + 1):(nsteps + 1)] <- 0
  cashflow[!exercise, t] <- 0
  
  # print(t)
}}

# Calculating option price by discounting all cashflows back to time 0
discounted_cashflows <- numeric(npaths)

for (i in 1:npaths) {
  for (t in 1:(nsteps + 1)) {
    if (cashflow[i, t] > 0) {  # If cashflow is non-zero (option exercised)
      discounted_cashflows[i] <- cashflow[i, t] * discount^(t-1)
      break  # Stop after the first non-zero cashflow (option exercise)
    }
  }
}

# Option price is the average of the first column of cashflows
cat("Bermudan Call Option Price:", mean(discounted_cashflows))
################################################################################




# Plot 1
gbm_plot_lines <- 50
plot(1:(nsteps + 1), S[1, ], type = "l", col = "blue",
     xlab = "Time Step", ylab = "Asset Price",
     main = "Simulated Asset Price Paths", ylim = range(S))
for (i in 2:gbm_plot_lines) {
  lines(1:(nsteps + 1), S[i, ], col = rainbow(gbm_plot_lines)[i])
}


# Plot 2 (Same plot as plotting stock price, )
plot(1:(nsteps + 1), colMeans(payoff), type = "l", col = "red",
     xlab = "Time Step", ylab = "Average Payoff",
     main = "Average Payoff Over Time")

# Plot 3
# Choose a time step (e.g., halfway through)
t <- floor(nsteps / 2)
in_the_money <- which(payoff[, t] > 0)
X <- S[in_the_money, t]
Y <- cashflow[in_the_money, t + 1] * discount
reg <- lm(Y ~ X + I(X^2))
plot(X, Y, col = "blue", pch = 16, xlab = "Stock Price", ylab = "Continuation Value",
     main = paste("Regression at Time Step", t))
curve(predict(reg, data.frame(X = x)), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Data", "Regression Fit"),
       col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))

# Plot 4
exercise_decision <- apply(cashflow, 1, function(row) which(row > 0)[1])
plot(1:npaths, (exercise_decision - 1), pch = 16, col = "purple",
     xlab = "Path Index", ylab = "First Exercise Time Step",
     main = "Exercise Decisions Across Paths")


# # Plot 5
# price_convergence <- numeric()
# for (i in seq(1000, npaths, by = 1000)) {
#   partial_cashflow <- cashflow[1:i, ]
#   price_convergence <- c(price_convergence, mean(partial_cashflow[, 1]) * discount^(-1))
# }
# plot(seq(1000, npaths, by = 1000), price_convergence, type = "l", col = "green",
#      xlab = "Number of Paths", ylab = "Option Price",
#      main = "Option Price Convergence")
# abline(h = mean(discounted_cashflows), col = "red", lty = 2)
# legend("topright", legend = c("Converging Price", "Final Estimate"),
#        col = c("green", "red"), lty = c(1, 2))
################################################################################




















