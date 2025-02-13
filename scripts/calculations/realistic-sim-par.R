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
registerDoParallel(cores = detectCores() - 2)

# Function to perform the simulation
simulate_bayes_factor <- function(mu, r, BF_crit, repetitions, trial_start = 2, trial_end = 0) {
  print(paste("Simulation starts for mu =", mu, "BF_crit =", BF_crit, "r =", r))
  start <- Sys.time()
  # Write the column headers
  results <- vector("list", length = repetitions)
  # Simulation
  for (i in 1:repetitions) {
    x <- rnorm(trial_start, mean = mu, sd = 1)
    BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
    stop_count <- trial_start
    # sym
    while (is.na(BF) || is.null(BF) || (BF > (1 / BF_crit) && BF < BF_crit &&
      # If a trial_end is specified, stop the simulation after that many trials
      if (trial_end > 0) {
        stop_count < trial_end
      } else {
        TRUE
      })) {
      x <- c(x, rnorm(1, mean = mu, sd = 1))
      BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
      stop_count <- stop_count + 1
    }
    decision <- ifelse(BF < BF_crit, 0, 1)
    results[[i]] <- list(decision = decision, stop_count = stop_count, mu = mu, bf_crit = BF_crit, r = r, trial_start = trial_start, trial_end = trial_end)
  }
  end <- Sys.time()
  print(paste("Simulation done for mu =", mu, "BF_crit =", BF_crit, "r =", r, "trial_start =", trial_start, "trial_end =", trial_end, "in", end - start))
  return(rbindlist(results))
}

# Initialize variables
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 200
# Connect to database
con <- dbConnect(duckdb(), "data/results.duckdb")
dbExecute(con, "DROP TABLE IF EXISTS cauchy_sym")
dbExecute(con, "CREATE TABLE IF NOT EXISTS cauchy_sym (decision INTEGER, stop_count INTEGER, mu DOUBLE, bf_crit DOUBLE, r DOUBLE, trial_start INTEGER, trial_end INTEGER)")
# Perform parallel simulations between [2, INF] trials
foreach(mu = big_sim_mus) %dopar% {
  foreach(r = r_vals) %do% {
    foreach(BF_crit = BF_crits) %do% {
      simulate_bayes_factor(mu, r, BF_crit, repetitions)
    } -> results
    results <- rbindlist(results)
    dbWriteTable(con, "cauchy_sym", results, append = TRUE)
  }
}
print("Done with the parallel symmetrical Optional Stopping simulation between [2, INF]")
# Perform parallel simulations between [20,200] trials
foreach(mu = big_sim_mus) %dopar% {
  foreach(r = r_vals[2]) %do% {
    foreach(BF_crit = BF_crits[1]) %do% {
      simulate_bayes_factor(mu, r, BF_crit, repetitions, trial_start = 20, trial_end = 200)
    } -> results
    results <- rbindlist(results)
    dbWriteTable(con, "cauchy_sym", results, append = TRUE)
  }
}
dbDisconnect(con)
# # Stop the parallel backend
stopImplicitCluster()
print("Done with the parallel symmetrical Optional Stopping simulation between [20, 200]")
