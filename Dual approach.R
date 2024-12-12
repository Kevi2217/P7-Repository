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
set.seed(1234567890)
S0 <- 100 # Underlying asset price
r <- 0.06 # risk-free rate
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

# Backward induction to compute continuation values
for (t in (nsteps + 1):1) {
  if (!(t %in% exercise_points)) {
    cashflow[, t] <- payoff[, t] * discount^(t-1)
    next
  }
  if (t == (nsteps + 1)) {
    cashflow[, t] <- payoff[, t]
    next
  }
  
  in_the_money <- which(payoff[, t] > 0)
  continuation_value <- numeric(npaths)
  if (length(in_the_money) > 0) {
    X <- S[in_the_money, t]
    Y <- payoff[in_the_money, t + 1] * discount^(t)
    continuation_value[in_the_money] <- predict(lm(Y ~ I(X) + I(X^2)),
                                                newdata = data.frame(X = S[in_the_money, t]))
    exercise <- payoff[, t] >= continuation_value
    cashflow[exercise, t] <- payoff[exercise, t]
    cashflow[exercise, (t + 1):(nsteps + 1)] <- 0
    cashflow[!exercise, t] <- 0
  }
}

cashflow[is.na(cashflow)] = 0

# Calculating option price by discounting all cashflows back to time 0
discounted_cashflows <- numeric(npaths)
for (i in 1:npaths) {
  for (t in 1:(nsteps + 1)) {
    if (cashflow[i, t] > 0) {
      discounted_cashflows[i] <- cashflow[i, t] * discount^(t-1)
      break
    }
  }
}

primal_estimate <- mean(discounted_cashflows)
cat("Bermudan Call Option Price (Primal):", primal_estimate, "\n")

############################# Dual Approach ####################################
# Function to compute dual price
compute_dual_price <- function(paths, exercise_values, continuation_values, payoff_matrix, discount, dt) {
  dual_values <- numeric(npaths)
  
  for (i in 1:paths) {
    martingale <- 0
    dual_value <- payoff_matrix[i, nsteps + 1] * discount^(nsteps)
    
    for (t in (nsteps):1) {
      in_the_money <- (payoff_matrix[i, t] > 0)
      if (in_the_money) {
        martingale <- martingale * discount + payoff_matrix[i, t]
        dual_value <- max(dual_value, martingale + payoff_matrix[i, t] - continuation_values[i, t])
      }
    }
    
    dual_values[i] <- dual_value
  }
  
  return(mean(dual_values))
}

# Compute continuation values from Longstaff-Schwartz
continuation_matrix <- matrix(0, nrow = npaths, ncol = nsteps + 1)
for (t in 1:nsteps) {
  in_the_money <- payoff[, t] > 0
  if (sum(in_the_money) > 0) {
    X <- S[in_the_money, t]
    Y <- payoff[in_the_money, t + 1] * discount^(t)
    continuation_matrix[in_the_money, t] <- predict(lm(Y ~ I(X) + I(X^2)),
                                                    newdata = data.frame(X = S[in_the_money, t]))
  }
}

# Compute dual price
dual_estimate <- compute_dual_price(npaths, cashflow, continuation_matrix, payoff, discount, dt)
cat("Bermudan Call Option Price (Dual):", dual_estimate, "\n")

################################################################################

# Comparing primal and dual estimates
cat("Primal Estimate:", primal_estimate, "\n")
cat("Dual Estimate:", dual_estimate, "\n")
