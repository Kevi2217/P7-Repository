
price_conv_function <- function (S0, r, nsteps, T, dt, K, sigma, discount, convergence_step) {
set.seed(1234567890)
price_convergence_list <- numeric(length(convergence_step))
avg_optinal_stopping <- numeric(length(convergence_step))
count_optimal_stopping <- numeric(length(convergence_step))
index <- 1
for (N in convergence_step) {
  # Simulating stock price
  dW <- matrix(rnorm(N * nsteps, mean = 0, sd = sqrt(dt)), ncol = nsteps)
  S <- matrix(0, nrow = N, ncol = nsteps + 1)
  
  S[, 1] <- S0
  for (t in 1:nsteps) {
    S[, t + 1] <- S[, t] * exp((r - 0.5 * sigma^2) * dt + sigma * dW[, t])
  }
  
  # Payoff matrix for a Bermudan Call option
  payoff <- pmax(S - K, 0)
  # Initialize cashflow matrix
  cashflow <- matrix(nrow = N, ncol = (nsteps + 1))
  
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
    continuation_value <- numeric(N)
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
    # print(t)
  }
  
  cashflow[is.na(cashflow)] = 0
  
  # Calculating option price by discounting all cashflows back to time 0
  discounted_cashflows <- numeric(N)
  
  for (i in 1:N) {
    for (t in 1:(nsteps + 1)) {
      if (cashflow[i, t] > 0) {  # If cashflow is non-zero (option exercised)
        discounted_cashflows[i] <- cashflow[i, t] * discount^(t-1)
        break  # Stop after the first non-zero cashflow (option exercise)
      }
    }
  }
price_convergence_list[index] <- mean(discounted_cashflows)
avg_optinal_stopping[index] <- mean((apply(cashflow, 1, function(row) which(row > 0)[1]) - 1), na.rm = TRUE)
count_optimal_stopping[index] <- sum(!is.na(apply(cashflow, 1, function(row) which(row > 0)[1])))
index <- index + 1   
}
return(list(price_convergence_list, avg_optinal_stopping, count_optimal_stopping))
}