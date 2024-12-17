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
r <- 0.1 # risk-free rate
nsteps <- 12 # Time steps
npaths <- 500000 # Number of paths
T <- 1 # Maturity (years)
dt <- T/nsteps # Discritization (time)
K <- 110 # Strike price
sigma <- 0.1 # Volatility
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
# To keep track of last exercises discounted
disc_cashflow <- payoff[, nsteps + 1] * discount
# Initialize cashflow matrix
cashflow <- matrix(nrow = npaths, ncol = (nsteps + 1))

############################# Longstaff Approach ###############################
# Backward induction to compute continuation values (ORG)
# t follows col numbers NOT time. So for discounts t is changed to match the time
# of the column. (There are 13 columns and 12 times. Col 1 represents time 0).
for (t in (nsteps + 1):1) {
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
    Y <- disc_cashflow[in_the_money]

    # Compute Laguerre basis functions
    L0 <- exp(-X/2)
    L1 <- exp(-X/2)*(1-X)
    L2 <- exp(-X/2)*(1-2*X+X^2/2)
    L3 <- exp(-X/2)*(1-3*X+(3*X^2/2)-(X^3/6))
    continuation_value[in_the_money] <- predict(lm(Y ~ L0 + L1 + L2 + L3), newdata = data.frame(L0 = exp(-X/2),
                                                                      L1 = exp(-X/2)*(1-X),
                                                                      L2 = exp(-X/2)*(1-2*X+X^2/2)),
                                                                      L3 = exp(-X/2)*(1-3*X+(3*X^2/2)-(X^3/6)))
    # Exercise or continue decision
    exercise <- payoff[, t] > continuation_value
    cashflow[exercise, t] <- payoff[exercise, t]
    cashflow[exercise, (t + 1):(nsteps + 1)] <- 0
    cashflow[!exercise, t] <- 0
    disc_cashflow <- ifelse(exercise, payoff[, t],
                                 disc_cashflow)
  }
    if (t != 1) {
      disc_cashflow <- disc_cashflow * discount
    }
  }

cashflow[is.na(cashflow)] = 0

# Option price is the average of the first column of cashflows
# exercise_matrix <- ifelse(cashflow != 0, 1, 0) %>%
#   colSums()

cat("Bermudan Call Option Price:", mean(disc_cashflow),
    (round(sum(!is.na(apply(cashflow, 1, function(row) which(row > 0)[1]))) / npaths * 100, 2)), "%")

######################### Black-Scholes Formula ################################
# Define the Black-Scholes formula for European Call and Put options
black_scholes <- function(S, K, r, T, sigma, type = "call") {
  d1 <- (log(S / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  if (type == "call") {
    price <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else if (type == "put") {
    price <- K * exp(-r * T) * pnorm(-d2) - S * pnorm(-d1)
  } else {
    stop("Invalid option type. Choose 'call' or 'put'.")
  }
  return(price)
}

# Calculate European Call and Put prices
call_price <- black_scholes(S = S0, K = K, r = r, T = T, sigma = sigma, type = "call")
put_price <- black_scholes(S = S0, K = K, r = r, T = T, sigma = sigma, type = "put")

# Output results
cat("European Call Option Price (Black-Scholes):", call_price, "\n")
cat("European Put Option Price (Black-Scholes):", put_price, "\n")

# Plot 1
S_plot <- as.data.frame(S[70020:7000, ]) %>% 
  mutate(observation = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "time", 
               values_to = "value") %>%
  mutate(time = as.numeric(factor(time)))

S_plot <- na.omit(S_plot)

gbm_plot_lines <- 20
plot(1:(nsteps + 1), S[1, ], type = "l", col = "blue",
     xlab = "Time Step", ylab = "Asset Price",
     main = "Simulated Asset Price Paths", ylim = range(S))
for (i in 2:gbm_plot_lines) {
  lines(1:(nsteps + 1), S[i, ], col = rainbow(gbm_plot_lines)[i])
}

# data.frame(time = unlist(apply(cashflow, 1, function(row) {which(row > 0)})) - 1)
# Plot 2
ggplot(data.frame(time = unlist(apply(cashflow, 1, function(row) {which(row > 0)})) - 1),
       aes(x = time)) +
  geom_bar(
    fill = "lightblue", 
    color = "black",  # Outline color
    size = 0.2,  # Thin outline
    stat = "count",  # Count the occurrences
    boundary = 0.5
  ) +
  scale_x_continuous(
    breaks = 0:12,  # Specify ticks from 1 to 12
    limits = c(0, 13)  # Ensure the axis shows only the range 1 to 12
  ) +
  labs(
    title = "Histogram of Stopping Times",
    x = "Time",
    y = "Frequency"
  )

# Plot 3
source("price_convergence_func.R")
convergence_step <- seq(from = 10000, to = 1000000, by = 10000)
price_convergence <- price_conv_function(S0, r, nsteps, T, dt,
                                         K, sigma, discount, convergence_step)
price_sd <- price_convergence[[4]] / sqrt(convergence_step)

cap <- 5000
ggplot(as.data.frame(cbind(convergence_step, price_convergence[[1]])),
       aes(x = convergence_step, y = price_convergence[[1]])) +
  geom_line(color = "blue") +
  geom_segment(aes(x = convergence_step, xend = convergence_step, 
                   y = price_convergence[[1]] - price_sd, 
                   yend = price_convergence[[1]] + price_sd),
               size = 0.2) +  # Thickness of the line
  # Horizontal caps at the top
  geom_segment(aes(x = convergence_step - cap,
                   xend = convergence_step + cap,
                   y = price_convergence[[1]] + price_sd, 
                   yend = price_convergence[[1]] + price_sd),
               size = 0.2) +  # Thickness of the cap line
  # Horizontal caps at the bottom
  geom_segment(aes(x = convergence_step - cap,
                   xend = convergence_step + cap,
                   y = price_convergence[[1]] - price_sd, 
                   yend = price_convergence[[1]] - price_sd),
               size = 0.2) +
  labs(
    x = "Number of Paths",
    y = "Option Price",
    title = "Option Price Based on Number of Paths"
  )

# Avg Exercise time over nPaths but without SE
# ggplot(as.data.frame(cbind(convergence_step, price_convergence[[2]])),
#        aes(x = convergence_step, y = price_convergence[[2]])) +
#   geom_line(color = "blue") +
#   labs(
#     x = "Number of Paths",
#     y = "Average Exercise Time",
#     title = "Average Exercise Time Based on Number of Paths"
#   )

# Plot 4
source("strike_price.R")
K_price <- seq(from = 0, to = 150, by = 5) 
K_conv <- strike_price_function(S0, r, nsteps, npaths, T,
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
# K_conv[[3]]/npaths * 100

# Plot 5
source("sigma_func.R")
sigma_steps <- seq(from = 0, to = 1, by = 0.05)
sigma_conv <- sigma_function(S0, r, nsteps, npaths, T,
                             dt, K, discount, sigma_steps)

ggplot(as.data.frame(cbind(sigma_steps, sigma_conv[[1]])),
       aes(x = sigma_steps, y = sigma_conv[[1]])) +
  geom_point(color = "blue") +
  labs(
    x = "Sigma",
    y = "Option Price",
    title = "Option Price Based on Volatility"
  ) 

ggplot(as.data.frame(cbind(sigma_steps, sigma_conv[[2]])),
       aes(x = sigma_steps, y = sigma_conv[[2]])) +
  geom_point(color = "blue") +
  labs(
    x = "Sigma",
    y = "Average Exercise Time",
    title = "Average Exercise Time Based on Volatility"
  )

# Plot 6
source("riskfree_func.R")
r_steps <- seq(from = 0, to = 0.5, by = 0.05)
r_conv <- riskfree_function(S0, nsteps, npaths, T, dt,
                            K, sigma, discount, r_steps)

ggplot(as.data.frame(cbind(r_steps, r_conv[[1]])),
       aes(x = r_steps, y = r_conv[[1]])) +
  geom_point(color = "blue") +
  labs(
    x = "Risk-Free Rate",
    y = "Option Price",
    title = "Option Price Based on Risk-Free Rate"
  ) 

ggplot(as.data.frame(cbind(r_steps, r_conv[[2]])),
       aes(x = r_steps, y = r_conv[[2]])) +
  geom_point(color = "blue") +
  labs(
    x = "Risk-Free Rate",
    y = "Average Exercise Time",
    title = "Average Exercise Time Based on Risk-Free Rate"
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
for (i in (nsteps+1):1) {
  for (j in 1:i) {
    if (i == (nsteps+1)) {
      stop_tree[, i] <- TRUE
      next
    }
    
    # Continuation value not discounted more than once, since it is just from the 1 time step discount to now
    continuation_value <- exp(-r * dt) * (p * payoff_tree[j + 1, i + 1] + (1 - p) * payoff_tree[j, i + 1])
    immediate_value <- payoff_tree[j, i]
    stop_tree[j, i] <- immediate_value >= continuation_value  # Record stopping decision
  }
}
stop_tree[payoff_tree <= 0] <- FALSE
# Initialize the option value tree at maturity
V_tree[, nsteps + 1] <- payoff_tree[, nsteps + 1]

# Backward induction for option price
for (i in nsteps:1) {
  for (j in 1:i) {
    continuation_value <- exp(-r * dt) * (p * V_tree[j + 1, i + 1] + (1 - p) * V_tree[j, i + 1])
    if (stop_tree[j, i]) {
      V_tree[j, i] <- payoff_tree[j, i]  # Immediate value if stopping
    } else {
      V_tree[j, i] <- continuation_value  # Continuation value otherwise
    }
  }
}

# Option price at time 0
option_price_tree <- V_tree[1, 1]
cat("Value of option (binomial model):", round(option_price_tree, 4), "\n")

