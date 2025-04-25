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
registerDoParallel(cores = detectCores() - 1)


#' Simulates Fixed Stopping from a sample size start and end.
simulate_fixed_size <- function(mu, r, BF_crit, repetitions, trial_start = 2, trial_end = 1000) {
    print(paste("Simulation starts with", repetitions, "repetitions", "for mu =", mu, "BF_crit =", BF_crit, "r =", r, "trial_start =", trial_start, "trial_end =", trial_end))
    start <- Sys.time()
    # Write the column headers
    results <- vector("list", length = (trial_end - trial_start + 1) * repetitions / 5)
    for (i in 1:repetitions) {
        # Simulation
        trial_count <- trial_start
        # fixed simulation
        while (trial_count <= trial_end) {
            x <- rnorm(trial_count, mean = mu, sd = 1)
            BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
            if (is.na(BF) || is.null(BF)) {
                next
            }
            # decide for 'H0', 'H1' or 'indecisive'
            decision <- if (BF < (1 / BF_crit)) {
                0
            } else if (BF > BF_crit) {
                1
            } else {
                2
            }
            results[[(i - 1) * trial_end + (trial_count - 1)]] <- list(decision = decision, trial_count = trial_count, bf = BF, mu = mu, bf_crit = BF_crit, r = r, trial_start = trial_start, trial_end = trial_end)
            trial_count <- trial_count + 5
        }
    }
    end <- Sys.time()
    print(paste("Simulation done with", repetitions, "repetitions", "for mu =", mu, "BF_crit =", BF_crit, "r =", r, "trial_start =", trial_start, "trial_end =", trial_end, "in", end - start))
    rbindlist(results)
}

# Initialize variables
big_sim_mus <- c(seq(0, 0.09, 0.01), seq(0.1, 1, 0.05))
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3)
repetitions <- 10000
chunk_size <- 500
# chunks need to be integers!
chunks <- 1:(repetitions / chunk_size)
# Connect to database
# con <- dbConnect(duckdb(), "data/hacking-bayes.duckdb") # WIP!!!!!!
# dbExecute(con, "DROP TABLE IF EXISTS cauchy_sym_fixed_size")
# dbExecute(con, "CREATE TABLE IF NOT EXISTS cauchy_sym_fixed_size (decision INTEGER, trial_count INTEGER, bf DOUBLE, mu DOUBLE, bf_crit DOUBLE, r DOUBLE, trial_start INTEGER, trial_end INTEGER)")
# # Perform parallel simulations between [2,1000] fixed sample sizes
# foreach(mu = big_sim_mus) %do% {
#     foreach(r = r_vals[2]) %do% {
#         foreach(BF_crit = BF_crits) %do% {
#             foreach(chunks) %dopar% {
#                 simulate_fixed_size(mu, r, BF_crit, chunk_size, trial_start = 2, trial_end = 1000)
#             } -> results
#             results <- rbindlist(results)
#             dbWriteTable(con, "cauchy_sym_fixed_size", results, append = TRUE)
#         }
#     }
# }
# dbDisconnect(con)

con <- dbConnect(duckdb(), "data/hacking-bayes.duckdb") # WIP!!!!!!
dbExecute(con, "DROP TABLE IF EXISTS cauchy_sym_fixed_size_r")
dbExecute(con, "CREATE TABLE IF NOT EXISTS cauchy_sym_fixed_size_r (decision INTEGER, trial_count INTEGER, bf DOUBLE, mu DOUBLE, bf_crit DOUBLE, r DOUBLE, trial_start INTEGER, trial_end INTEGER)")
repetitions <- 10000
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3)
chunk_size <- 500
# chunks need to be integers!
chunks <- 1:(repetitions / chunk_size)
foreach(mu = big_sim_mus) %do% {
    foreach(r = r_vals) %do% {
        foreach(BF_crit = BF_crits) %do% {
            foreach(chunks) %dopar% {
                simulate_fixed_size(mu, r, BF_crit, chunk_size, trial_start = 2, trial_end = 1000)
            } -> results
            results <- rbindlist(results)
            dbWriteTable(con, "cauchy_sym_fixed_size_r", results, append = TRUE)
        }
    }
}
dbDisconnect(con)

con <- dbConnect(duckdb(), "data/hacking-bayes.duckdb") # WIP!!!!!!
dbExecute(con, "DROP TABLE IF EXISTS cauchy_sym_fixed_size_bf_crit")
dbExecute(con, "CREATE TABLE IF NOT EXISTS cauchy_sym_fixed_size_bf_crit (decision INTEGER, trial_count INTEGER, bf DOUBLE, mu DOUBLE, bf_crit DOUBLE, r DOUBLE, trial_start INTEGER, trial_end INTEGER)")
repetitions <- 10000
BF_crits <- c(3, 6, 10)
chunk_size <- 500
# chunks need to be integers!
chunks <- 1:(repetitions / chunk_size)
foreach(mu = big_sim_mus) %do% {
    foreach(r = r_vals[2]) %do% {
        foreach(BF_crit = BF_crits) %do% {
            foreach(chunks) %dopar% {
                simulate_fixed_size(mu, r, BF_crit, chunk_size, trial_start = 2, trial_end = 1000)
            } -> results
            results <- rbindlist(results)
            dbWriteTable(con, "cauchy_sym_fixed_size_bf_crit", results, append = TRUE)
        }
    }
}
dbDisconnect(con)
# # Stop the parallel backend
stopImplicitCluster()
