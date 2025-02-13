# Set working directory
setwd("/home/linus/git/hacking-bayes")

# Clear workspace
rm(list = ls())
gc()

# Import necessary libraries
library(BayesFactor)
library(foreach)
library(doParallel)
library(data.table)
library(duckdb)

# Register parallel backend
registerDoParallel(cores = detectCores() - 1)

# Function to perform the simulation
simulate_bayes_factor <- function(mu, r, BF_crit, repetitions) {
  print(paste("Simulation starts for mu =", mu, "BF_crit =", BF_crit, "r =", r))
  start <- Sys.time()
  # Write the column headers
  results <- vector("list", length = repetitions)
  # Simulation
  for (i in 1:repetitions) {
    x <- rnorm(2, mean = mu, sd = 1)
    BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
    stop_count <- 1
    # sym
    while (is.na(BF) || is.null(BF) || (BF > (1 / BF_crit) && BF < BF_crit)) {
      x <- c(x, rnorm(1, mean = mu, sd = 1))
      BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
      stop_count <- stop_count + 1
    }
    decision <- ifelse(BF < BF_crit, 0, 1)
    results[[i]] <- list(decision = decision, stop_count = stop_count, mu = mu, bf_crit = BF_crit, r = r)
  }
  end <- Sys.time()
  print(paste("Simulation done for mu =", mu, "BF_crit =", BF_crit, "r =", r, "in", end - start))
  return(rbindlist(results))
}

# Initialize variables
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000
# Connect to database
con <- dbConnect(duckdb(), "data/results.duckdb")
dbExecute(con, "DROP TABLE IF EXISTS cauchy_sym")
dbExecute(con, "CREATE TABLE IF NOT EXISTS cauchy_sym (decision INTEGER, stop_count INTEGER, mu DOUBLE, bf_crit DOUBLE, r DOUBLE)")
# Perform parallel simulations
for (mu in big_sim_mus) {
  foreach(r = r_vals) %:%
    foreach(BF_crit = BF_crits) %dopar% {
      simulate_bayes_factor(mu, r, BF_crit, repetitions)
    } -> results
  results <- rbindlist(unlist(results, recursive = FALSE))
  dbWriteTable(con, "cauchy_sym", results, append = TRUE)
}
dbDisconnect(con)
# # Stop the parallel backend
stopImplicitCluster()
