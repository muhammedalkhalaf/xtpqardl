set.seed(42)
n_countries <- 10
n_years <- 30
years <- 1990:(1990 + n_years - 1)
beta1 <- 0.8; beta2 <- 0.5; gamma1 <- 0.4; gamma2 <- 0.3
pqardl_sample <- data.frame()
for (i in 1:n_countries) {
  mu_i <- rnorm(1, mean = 5, sd = 1)
  rho_i <- runif(1, min = -0.6, max = -0.2)
  y <- x1 <- x2 <- numeric(n_years)
  x1[1] <- rnorm(1, mean = 3, sd = 0.5)
  x2[1] <- rnorm(1, mean = 2, sd = 0.5)
  y[1] <- mu_i + beta1 * x1[1] + beta2 * x2[1] + rnorm(1, sd = 0.3)
  for (t in 2:n_years) {
    x1[t] <- x1[t-1] + rnorm(1, mean = 0.05, sd = 0.2)
    x2[t] <- x2[t-1] + rnorm(1, mean = 0.03, sd = 0.15)
    ect <- y[t-1] - mu_i - beta1 * x1[t-1] - beta2 * x2[t-1]
    d_x1 <- x1[t] - x1[t-1]; d_x2 <- x2[t] - x2[t-1]
    u <- rnorm(1, sd = 0.3 * (1 + 0.2 * abs(ect)))
    d_y <- rho_i * ect + gamma1 * d_x1 + gamma2 * d_x2 + u
    y[t] <- y[t-1] + d_y
  }
  cd <- data.frame(country = factor(paste0("Country_", LETTERS[i])), year = years, y = y, x1 = x1, x2 = x2)
  cd$L_y <- c(NA, y[-n_years]); cd$d_y <- c(NA, diff(y)); cd$d_x1 <- c(NA, diff(x1)); cd$d_x2 <- c(NA, diff(x2))
  pqardl_sample <- rbind(pqardl_sample, cd)
}
pqardl_sample <- pqardl_sample[!is.na(pqardl_sample$d_y), ]
rownames(pqardl_sample) <- NULL
save(pqardl_sample, file = "data/pqardl_sample.rda", compress = "xz")
cat("Data saved successfully\n")
