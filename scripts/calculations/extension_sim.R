#
# Combined fixed-N + optional-stopping simulation for large effect sizes (mu = 2, 3)
#
setwd(".")

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

# ---- Simulation functions ----------------------------------------------

#' Fixed-N sampling: draw exactly trial_count observations, compute BF once.
simulate_fixed_size <- function(mu, r, BF_crit, repetitions, trial_start = 2, trial_end = 50, high_density = TRUE) {
  print(paste("[fixed] mu =", mu, "r =", r, "BF_crit =", BF_crit, "reps =", repetitions, "N in [", trial_start, ",", trial_end, "]"))
  start <- Sys.time()
  results <- vector("list", (trial_end - trial_start + 1) * repetitions)
  idx <- 1
  for (i in 1:repetitions) {
    trial_count <- trial_start
    while (trial_count <= trial_end) {
      x <- rnorm(trial_count, mean = mu, sd = 1)
      BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
      if (is.na(BF) || is.null(BF)) {
        trial_count <- trial_count + 1
        next
      }
      decision <- if (BF < (1 / BF_crit)) 0 else if (BF > BF_crit) 1 else 2
      results[[idx]] <- list(decision = decision, trial_count = trial_count, bf = BF,
                              mu = mu, bf_crit = BF_crit, r = r,
                              trial_start = trial_start, trial_end = trial_end)
      idx <- idx + 1
      trial_count <- if (high_density && trial_count < 50) trial_count + 1 else trial_count + 5
    }
  }
  end <- Sys.time()
  print(paste("[fixed] done in", round(difftime(end, start, units = "secs"), 1), "sec"))
  rbindlist(results[1:(idx - 1)])
}

#' Optional stopping: keep sampling one at a time until BF crosses criterion or trial_end is hit.
simulate_optional_stopping <- function(mu, r, BF_crit, repetitions, trial_start = 2, trial_end = 200) {
  print(paste("[optional] mu =", mu, "r =", r, "BF_crit =", BF_crit, "reps =", repetitions, "N in [", trial_start, ",", trial_end, "]"))
  start <- Sys.time()
  results <- vector("list", repetitions)
  for (i in 1:repetitions) {
    x <- rnorm(trial_start, mean = mu, sd = 1)
    BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
    stop_count <- trial_start
    while ((is.na(BF) || is.null(BF) ||
            (BF > (1 / BF_crit) && BF < BF_crit)) && stop_count < trial_end) {
      x <- c(x, rnorm(1, mean = mu, sd = 1))
      BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
      stop_count <- stop_count + 1
    }
    decision <- if (is.na(BF) || is.null(BF)) 2 else if (BF < (1 / BF_crit)) 0 else if (BF > BF_crit) 1 else 2
    results[[i]] <- list(decision = decision, stop_count = stop_count, mu = mu,
                          bf_crit = BF_crit, r = r, trial_start = trial_start, trial_end = trial_end)
  }
  end <- Sys.time()
  print(paste("[optional] done in", round(difftime(end, start, units = "secs"), 1), "sec"))
  rbindlist(results)
}

# ---- Config --------------------------------------------------------------
mus      <- c(2, 3)
r_vals   <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3)

fixed_trial_start    <- 2
fixed_trial_end      <- 50
fixed_repetitions    <- 10000
fixed_chunk_size     <- 500
fixed_chunks         <- 1:(fixed_repetitions / fixed_chunk_size)

optional_trial_start <- 20
optional_trial_end   <- 200
optional_repetitions <- 20000
optional_chunk_size  <- 2000
optional_chunks      <- 1:(optional_repetitions / optional_chunk_size)

# ---- Run -------------------------------------------------------------

con <- dbConnect(duckdb(), "data/large_effects.duckdb")

dbExecute(con, "DROP TABLE IF EXISTS fixed_size_large_effect")
dbExecute(con, "CREATE TABLE fixed_size_large_effect (decision INTEGER, trial_count INTEGER, bf DOUBLE, mu DOUBLE, bf_crit DOUBLE, r DOUBLE, trial_start INTEGER, trial_end INTEGER)")

dbExecute(con, "DROP TABLE IF EXISTS optional_stopping_large_effect")
dbExecute(con, "CREATE TABLE optional_stopping_large_effect (decision INTEGER, stop_count INTEGER, mu DOUBLE, bf_crit DOUBLE, r DOUBLE, trial_start INTEGER, trial_end INTEGER)")

foreach(mu = mus) %do% {
  foreach(r = r_vals) %do% {
    foreach(BF_crit = BF_crits) %do% {
      # Fixed-N
      foreach(fixed_chunks) %dopar% {
        simulate_fixed_size(mu, r, BF_crit, fixed_chunk_size, trial_start = fixed_trial_start, trial_end = fixed_trial_end)
      } -> fixed_results
      dbWriteTable(con, "fixed_size_large_effect", rbindlist(fixed_results), append = TRUE)

      # Optional stopping
      foreach(optional_chunks) %dopar% {
        simulate_optional_stopping(mu, r, BF_crit, optional_chunk_size, trial_start = optional_trial_start, trial_end = optional_trial_end)
      } -> optional_results
      dbWriteTable(con, "optional_stopping_large_effect", rbindlist(optional_results), append = TRUE)
    }
  }
}

dbDisconnect(con)
stopImplicitCluster()

print("Done: fixed-N and optional-stopping simulations for mu = 2, 3")