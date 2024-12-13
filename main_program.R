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
  # Only consider in-the-money paths (Continuation_value is defined previous to definition to ensure indices align with cashflow-matrix)
  in_the_money <- which(payoff[, t] > 0)
  continuation_value <- numeric(npaths)
  if (length(in_the_money) > 0) {
    
  # Fit regression to estimate continuation value
  X <- S[in_the_money, t]
  Y <- payoff[in_the_money, t + 1] * discount^(t)
  continuation_value[in_the_money] <- predict(lm(Y ~ I(X) + I(X^2)),
                                newdata = data.frame(X = S[in_the_money, t]))
  # Exercise or continue decision
  exercise <- payoff[, t] >= continuation_value
  cashflow[exercise, t] <- payoff[exercise, t]
  cashflow[exercise, (t + 1):(nsteps + 1)] <- 0
  cashflow[!exercise, t] <- 0
  }
print(t)
}

cashflow[is.na(cashflow)] = 0

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
# Extract 20 paths for plotting
plot_paths <- 20
S_plot <- S[1:plot_paths, ]

# Convert data to long format for ggplot2
time_steps <- 0:nsteps
stock_prices <- data.frame(time = rep(time_steps, each = plot_paths),
                           path = rep(1:plot_paths, times = length(time_steps)),
                           price = as.vector(S_plot))

# Plot stock prices
p <- ggplot(stock_prices, aes(x = time, y = price, group = path, color = factor(path))) +
  geom_line() +
  scale_x_continuous(breaks = time_steps) +
  labs(title = "Simulated Stock Price Paths",
       x = "Time",
       y = "Stock Price",
       color = "Path") +
  theme(legend.position="none")

print(p)

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


# Plot 5
source("price_convergence_func.R")
convergence_step <- seq(from = 10000, to = 1000000, by = 10000) # c(100, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
price_convergence <- price_conv_function(S0, r, nsteps, T, dt,
                                         K, sigma, discount, convergence_step)

ggplot(as.data.frame(cbind(convergence_step, price_convergence[[1]])),
       aes(x = convergence_step, y = price_convergence[[1]])) +
  geom_line(color = "blue") +
  labs(
    x = "Number of Paths",
    y = "Option Price",
    title = "Option Price Based on Number of Paths"
  ) 
ggplot(as.data.frame(cbind(convergence_step, price_convergence[[2]])),
       aes(x = convergence_step, y = price_convergence[[2]])) +
  geom_line(color = "blue") +
  labs(
    x = "Number of Paths",
    y = "Average Exercise Time",
    title = "Average Exercise Time Based on Number of Paths"
  )

# Plot 6
source("strike_price.R")
K_price <- seq(from = 0, to = 150, by = 5) 
K_conv <- strike_price_function(S0, r, nsteps, 500000, T,
                                dt, sigma, discount, K_price)

ggplot(as.data.frame(cbind(K_price, K_conv[[1]])),
       aes(x = K_price, y = K_conv[[1]])) +
  geom_point(color = "blue") +
  labs(
    x = "Strike Price",
    y = "Option Price",
    title = "Option Price Based on Strike Price"
  ) 

ggplot(as.data.frame(cbind(K_price, K_conv[[2]])),
       aes(x = K_price, y = K_conv[[2]])) +
  geom_point(color = "blue") +
  labs(
    x = "Strike Price",
    y = "Average Exercise Time",
    title = "Average Exercise Time Based on Strike Price"
  ) 
# K_conv[[3]]/500000 * 100

# Plot 7
source("sigma_func.R")
sigma_steps <- seq(from = 0, to = 1, by = 0.05)
sigma_conv <- sigma_function(S0, r, nsteps, 500000, T,
                             dt, K, discount, sigma_steps)

ggplot(as.data.frame(cbind(sigma_steps, sigma_conv[[1]])),
       aes(x = sigma_steps, y = sigma_conv[[1]])) +
  geom_point(color = "blue") +
  labs(
    x = "Sigma",
    y = "Option Price",
    title = "Option Price Based on Sigma"
  ) 

ggplot(as.data.frame(cbind(sigma_steps, sigma_conv[[2]])),
       aes(x = sigma_steps, y = sigma_conv[[2]])) +
  geom_point(color = "blue") +
  labs(
    x = "Sigma",
    y = "Average Exercise Time",
    title = "Average Exercise Time Based on Sigma"
)

# Plot 8
source("riskfree_func.R")
r_steps <- seq(from = 0, to = 0.3, by = 0.02)
riskfree_conv <- riskfree_function(S0, nsteps, 500000, T, dt,
                                   K, sigma, discount, r_steps)

ggplot(as.data.frame(cbind(r_steps, riskfree_conv[[1]])),
       aes(x = r_steps, y = riskfree_conv[[1]])) +
  geom_point(color = "blue") +
  labs(
    x = "Riskfree-Rate",
    y = "Option Price",
    title = "Option Price Based on Riskfree-Rate"
  ) 

ggplot(as.data.frame(cbind(r_steps, riskfree_conv[[2]])),
       aes(x = r_steps, y = riskfree_conv[[2]])) +
  geom_point(color = "blue") +
  labs(
    x = "Riskfree-Rate",
    y = "Average Exercise Time",
    title = "Average Exercise Time Based on Riskfree-Rate"
  )



########################### Binomial Model ###############################
# Calculate up and down factors
u <- exp(sigma * sqrt(dt))        # Up factor
d <- 1 / u                        # Down factor
p <- (exp(r * dt) - d) / (u - d)  # Risk-neutral probability

# Initialize matrices
S_tree <- matrix(0, nrow = nsteps + 1, ncol = nsteps + 1)  # Stock price tree
V_tree <- matrix(0, nrow = nsteps + 1, ncol = nsteps + 1)  # Option value tree
stop_tree <- matrix(FALSE, nrow = nsteps + 1, ncol = nsteps + 1)  # Stopping decision tree

# Build the stock price tree
for (i in 0:nsteps) {
  for (j in 0:i) {
    S_tree[j + 1, i + 1] <- S0 * (u^j) * (d^(i - j))
  }
}

# payoff tree (intrinsic value of exercising at each node)
payoff_tree <- pmax(S_tree - K, 0)

# Backward induction
for (i in (nsteps):1) {
  for (j in 1:i) {
    # Continuation value not discounted more than once, since it is just from the 1 time step discount to now
    continuation_value <- exp(-r * dt) * (p * payoff_tree[j + 1, i + 1] + (1 - p) * payoff_tree[j, i + 1])
    immediate_value <- payoff_tree[j, i]
    stop_tree[j, i] <- immediate_value >= continuation_value  # Record stopping decision
    # cat("i:", i," j:", j, "\n")
  }
}


# Result: Option value at the root
value_binomial <- payoff_tree[1, 1] # Forkert
cat("Fair value of the option (binomial model):", round(value_binomial, 4), "\n")

stop_df <- reshape2::melt(stop_tree, varnames = c("Step", "Node"), value.name = "Stop")
ggplot(stop_df, aes(x = Node, y = Step, fill = Stop)) +
  geom_tile() +
  labs(title = "Stopping Decisions in Binomial Tree")

















