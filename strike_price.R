
strike_price_function <- function (S0, r, nsteps, npaths, T, dt, sigma, discount, K_price) {
  set.seed(1234567890)
  return_list <- numeric(length(K_price))
  avg_optimal_stopping <- numeric(length(K_price))
  count_optimal_stopping <- numeric(length(K_price))
  index <- 1
  
  for (K in K_price) {
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
    
    return_list[index] <- mean(disc_cashflow)
    avg_optimal_stopping[index] <- mean((apply(cashflow, 1, function(row) which(row > 0)[1]) - 1), na.rm = TRUE)
    count_optimal_stopping[index] <- sum(!is.na(apply(cashflow, 1, function(row) which(row > 0)[1])))
    index <- index + 1   
  }
  return(list(return_list, avg_optimal_stopping, count_optimal_stopping))
}