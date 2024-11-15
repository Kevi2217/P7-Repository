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


####################################### GBM ####################################
# set.seed(21)
mu <- 1
n <- 50
dt <- 0.1
x0 <- 100


sigma <- seq(0.8, 1.8, by = 0.2)

# Random paths
x <- sapply(sigma, function(s) {
  exp((mu - s^2 / 2) * dt + s * rnorm(n, mean = 0, sd = sqrt(dt)))
})
x <- rbind(rep(1, length(sigma)), x)
x <- apply(x, 2, cumprod) * x0

# Convert data to long format
x_df <- as.data.frame(x)
x_df$time <- 0:n
x_long <- reshape2::melt(x_df, id.vars = "time", variable.name = "sigma", value.name = "x")

ggplot(x_long, aes(x = time, y = x, color = factor(sigma))) +
  geom_line() +
  scale_color_discrete(name = "Sigma", labels = round(sigma, 2)) +
  labs(
    x = expression(t),
    y = expression(x),
    title = expression("Geometric Brownian Motion with different variances" ~ (mu == 1))
  ) +
  theme_minimal()
################################################################################





















