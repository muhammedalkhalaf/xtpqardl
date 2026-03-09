# Generate simulated panel data for PQARDL package
# This script creates the pqardl_sample dataset

set.seed(42)

# Parameters
n_countries <- 10
n_years <- 30
years <- 1990:(1990 + n_years - 1)

# True parameters
beta1 <- 0.8    # Long-run effect of x1
beta2 <- 0.5    # Long-run effect of x2
gamma1 <- 0.4   # Short-run effect of d_x1
gamma2 <- 0.3   # Short-run effect of d_x2

# Initialize data
pqardl_sample <- data.frame()

for (i in 1:n_countries) {
  # Country-specific parameters
  mu_i <- rnorm(1, mean = 5, sd = 1)         # Fixed effect
  rho_i <- runif(1, min = -0.6, max = -0.2)  # ECT speed of adjustment
  
  # Initialize variables
  y <- numeric(n_years)
  x1 <- numeric(n_years)
  x2 <- numeric(n_years)
  
  # Initial values
  x1[1] <- rnorm(1, mean = 3, sd = 0.5)
  x2[1] <- rnorm(1, mean = 2, sd = 0.5)
  y[1] <- mu_i + beta1 * x1[1] + beta2 * x2[1] + rnorm(1, sd = 0.3)
  
  for (t in 2:n_years) {
    # Generate x variables (random walk with drift)
    x1[t] <- x1[t-1] + rnorm(1, mean = 0.05, sd = 0.2)
    x2[t] <- x2[t-1] + rnorm(1, mean = 0.03, sd = 0.15)
    
    # Error correction term
    ect <- y[t-1] - mu_i - beta1 * x1[t-1] - beta2 * x2[t-1]
    
    # Changes in x
    d_x1 <- x1[t] - x1[t-1]
    d_x2 <- x2[t] - x2[t-1]
    
    # Error term with some heteroskedasticity based on quantile
    u <- rnorm(1, sd = 0.3 * (1 + 0.2 * abs(ect)))
    
    # ECM equation
    d_y <- rho_i * ect + gamma1 * d_x1 + gamma2 * d_x2 + u
    
    y[t] <- y[t-1] + d_y
  }
  
  # Create data frame for this country
  country_data <- data.frame(
    country = factor(paste0("Country_", LETTERS[i])),
    year = years,
    y = y,
    x1 = x1,
    x2 = x2
  )
  
  # Add lagged and differenced variables
  country_data$L_y <- c(NA, y[-n_years])
  country_data$d_y <- c(NA, diff(y))
  country_data$d_x1 <- c(NA, diff(x1))
  country_data$d_x2 <- c(NA, diff(x2))
  
  pqardl_sample <- rbind(pqardl_sample, country_data)
}

# Remove first observation per country (NA due to differencing)
pqardl_sample <- pqardl_sample[!is.na(pqardl_sample$d_y), ]

# Reset row names
rownames(pqardl_sample) <- NULL

# Save to data directory
usethis::use_data(pqardl_sample, overwrite = TRUE)
