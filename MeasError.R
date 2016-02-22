library(dplyr)
library(leiv)

# Function to bootstrap sample covariance matrix S

boot_S <- function(x, y, iter = 999) 
{
  require(dplyr)
  data <- cbind(x, y)
  
  boot_fun <- function(data)
  {
    N <- nrow(data)
    samps <- sample(1:N, N, replace = TRUE)
    data_new <- data[samps, ]
    cov_est <- cov(data_new)
    data.frame(S11 = cov_est[1,1], S12 = cov_est[1,2], 
               S21 = cov_est[2,1], S22 = cov_est[2,2])
  }
  
  boot_reps <- replicate(iter, boot_fun(data), simplify = FALSE)
  boot_df <- bind_rows(boot_reps)
  boot_df
}

# Perform asymptotic Bayesian, Geometric Mean, OLS Bisector, and Orthogonal
# Regression.

all_reg <- function(x, y)
{
  require(dplyr)
  require(leiv)
  
  boot_df <- boot_S(x, y)
  quant_df <- boot_df %>%
    mutate(b1 = S12 / S11) %>%
    mutate(b2 = S22 / S12) %>%
    mutate(B = 0.5 * (b2 - 1 / b1)) %>%
    mutate(geom_est = sign(S12) * sqrt(b1 * b2)) %>%
    mutate(bis_est = tan(0.5 * (atan(b1) + atan(b2)))) %>%
    mutate(orth_est = B + sign(S12) * sqrt(B^2 + 1))
  
  data <- cbind(x, y)
  cov_est <- cov(data)
  S11 <- cov_est[1,1]; S12 <- cov_est[1,2]; 
  S21 <- cov_est[2,1]; S22 <- cov_est[2,2];
  b1 <- S12 / S11; b2 <- S22 / S12
  theta1 <- atan(b1); theta2 <- atan(b2)
  B <- 0.5 * (b2 - 1 / b1)
  
  leiv_fit <- leiv(y ~ x, probIntCalc = TRUE, level = 0.90, abs.tol = 1e-10)
  bayes_est <- leiv_fit@slope
  bayes_ci <- leiv_fit@slopeInt
  
  geom_est = sign(S12) * sqrt(b1 * b2)
  bis_est = tan(0.5 * (atan(b1) + atan(b2)))
  orth_est = B + sign(S12) * sqrt(B^2 + 1)
  
  geom_ci <- quantile(quant_df$geom_est, probs = c(0.05, 0.95))
  bis_ci <- quantile(quant_df$bis_est, probs = c(0.05, 0.95))
  orth_ci <- quantile(quant_df$orth_est, probs = c(0.05, 0.95))
  
  out <- data.frame(method = c('Bayes', 'Geometric', 'OLSBis', 'OrthReg'),
                    estimate = c(bayes_est, geom_est, bis_est, orth_est),
                    lower = c(bayes_ci[1], geom_ci[1], bis_ci[1], orth_ci[1]),
                    upper = c(bayes_ci[2], geom_ci[2], bis_ci[2], orth_ci[2]))
  out
}

gen_data <- function(n, sigma1, sigma2)
{
  xi_1 <- rnorm(n, 0, 1)
  x <- xi_1 + rnorm(n, 0, sigma1)
  y <- x + rnorm(n, 0, sigma2)
  data.frame(x = x, y = y)
}

n <- 20; sigma1 <- 0.50; sigma2 <- 0.10

reg_sim <- function(n, sigma1, sigma2)
{
  data <- gen_data(n, sigma1, sigma2)
  out_df <- all_reg(data$x, data$y)
  out_df
}

# Small simulation

res_list <- replicate(100, reg_sim(n, sigma1, sigma2), simplify = FALSE)

res_df <- bind_rows(res_list)

cover <- res_df %>%
  group_by(method) %>%
  summarise(cover = mean(lower < 1 & upper > 1))

# Try a different construction for the CI of the orthogonal estimate
# Based on Jackson and Dunleavy 1988

all_reg2 <- function(x, y)
{
  require(dplyr)
  require(leiv)
  
  data <- cbind(x, y)
  n <- nrow(data)
  cov_est <- cov(data)
  eig_vals <- eigen(cov_est)$values
  l1 <- eig_vals[1]
  l2 <- eig_vals[2]
  phi <- asin( (4 * pf(0.05, 1, n-2) / 
                  ((n - 2) * (l1 / l2 + l2 / l1 - 2)))^(1/2) ) / 2
  
  S11 <- cov_est[1,1]; S12 <- cov_est[1,2]; 
  S21 <- cov_est[2,1]; S22 <- cov_est[2,2];
  b1 <- S12 / S11; b2 <- S22 / S12
  theta1 <- atan(b1); theta2 <- atan(b2)
  B <- 0.5 * (b2 - 1 / b1)
  
  leiv_fit <- leiv(y ~ x, probIntCalc = TRUE, level = 0.90, abs.tol = 1e-10)
  bayes_est <- leiv_fit@slope
  bayes_ci <- leiv_fit@slopeInt
  
  orth_est = B + sign(S12) * sqrt(B^2 + 1)
  theta_hat <- atan(orth_est)
  orth_ci <- c(tan(theta_hat - phi), tan(theta_hat + phi))
  
  out <- data.frame(method = c('Bayes', 'OrthReg'),
                    estimate = c(bayes_est, orth_est),
                    lower = c(bayes_ci[1], orth_ci[1]),
                    upper = c(bayes_ci[2], orth_ci[2]))
  out
}


reg_sim2 <- function(n, sigma1, sigma2)
{
  data <- gen_data(n, sigma1, sigma2)
  out_df <- all_reg2(data$x, data$y)
  out_df
}

res_list <- replicate(100, reg_sim2(n, sigma1, sigma2), simplify = FALSE)

res_df <- bind_rows(res_list)

cover <- res_df %>%
  group_by(method) %>%
  summarise(cover = mean(lower < 1 & upper > 1))
