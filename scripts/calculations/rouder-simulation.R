setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# import bf
source("bayes-factor-functions.R")

# load libraries
library(data.table)
library(duckdb)

# fixed_mu <- 0.1
# mus <- c(0.1, 0.2, 0.5, 0.8, 1)
big_sim_mus <- seq(0, 1, 0.01)
fixed_sigma <- 1
sigmas <- c(0.01, 0.1, 1, 5, 10)
repetitions <- 20000
n_starts <- 1:15

# Create a table to store the results
con <- dbConnect(duckdb(), dbdir = "data/results.duckdb")
dbExecute(con, "DROP TABLE IF EXISTS bf_decision_threshold")
dbExecute(con, "CREATE TABLE bf_decision_threshold (decision INTEGER, stop_count INTEGER, mu DOUBLE, trial_start INTEGER, symmetrical BOOLEAN)")
dbDisconnect(con)
######################################################################
#### Simulation for asymmetrical Bayes Factor decision thresholds ####
####                                                              ####
######################################################################
for (k in n_starts) {
  results <- vector("list", length = length(big_sim_mus) * 20000)
  mu_counter <- 0
  for (mu in big_sim_mus) {
    # reset seed
    set.seed(NULL)
    # reset h0/h1_counts
    h0_count <- 0
    indecisive_count <- 0
    # run the simulation 20000 times
    for (i in 1:repetitions) {
      h1 <- rnorm(k, mean = mu, sd = fixed_sigma)
      bf <- bf01(bf10(length(h1), mean(h1), var(h1)))
      if (is.na(bf)) {
        bf <- 1.33
      }
      # run the simulation until the bayes factor is between 1/3 and 3
      stop_count <- 1
      while (bf < 3 && stop_count < 250) {
        h1 <- c(h1, rnorm(1, mean = mu, sd = 1))
        bf <- bf01(bf10(length(h1), mean(h1), var(h1)))
        stop_count <- stop_count + 1
      }
      # count the number of times the null hypothesis is true
      # and save the results for plotting
      if (bf > 3) {
        h0_count <- h0_count + 1
        results[[mu_counter * repetitions +
          h0_count + indecisive_count]] <- list(0, stop_count, mu, k, FALSE)
      } else {
        indecisive_count <- indecisive_count + 1
        results[[mu_counter * repetitions +
          h0_count + indecisive_count]] <- list(1, stop_count, mu, k, FALSE)
      }
    }
    mu_counter <- mu_counter + 1
    # print the results
    print(paste("Asymetrical simulation is done for mu =", mu))
    print(paste("H0: ", h0_count))
    print(paste("H1 or Indecisive: ", indecisive_count))
  }
  # Export results
  results <- rbindlist(results)
  setnames(results, c("decision", "stop_count", "mu", "trial_start", "symmetrical"))
  con <- dbConnect(duckdb(), dbdir = "data/results.duckdb")
  dbWriteTable(con, "bf_decision_threshold", results, append = TRUE)
}
#####################################################################
#### Simulation for symmetrical Bayes Factor decision thresholds ####
#####################################################################

for (k in n_starts) {
  results <- vector("list", length = length(big_sim_mus) * 20000)
  mu_counter <- 0
  for (mu in big_sim_mus) {
    # reset seed
    set.seed(NULL)
    # reset h0/h1_counts
    h0_count <- 0
    h1_count <- 0
    # run the simulation 20000 times
    for (i in 1:repetitions) {
      h1 <- rnorm(k, mean = mu, sd = fixed_sigma)
      bf <- bf01(bf10(length(h1), mean(h1), var(h1)))
      if (is.na(bf)) {
        bf <- 1.33
      }
      # run the simulation until the bayes factor is between 1/3 and 3
      stop_count <- 1
      while (bf < 3 && bf > 1 / 3) {
        h1 <- c(h1, rnorm(1, mean = mu, sd = 1))
        bf <- bf01(bf10(length(h1), mean(h1), var(h1)))
        stop_count <- stop_count + 1
      }
      # count the number of times the null hypothesis is true
      # and save the results for plotting
      if (bf > 3) {
        h0_count <- h0_count + 1
        results[[mu_counter * repetitions +
          h0_count + h1_count]] <- list(0, stop_count, mu, k, TRUE)
      } else {
        h1_count <- h1_count + 1
        results[[mu_counter * repetitions +
          h0_count + h1_count]] <- list(1, stop_count, mu, k, TRUE)
      }
      # print the results
      # print(paste("Decision: ", bf, "Stop count: ", stop_count))
    }
    mu_counter <- mu_counter + 1
    # print the results
    print(paste("Symmetrical simulation is done for mu =", mu))
    print(paste("H0: ", h0_count))
    print(paste("H1: ", h1_count))
  }
  # Export results
  results <- rbindlist(results)
  setnames(results, c("decision", "stop_count", "mu", "trial_start", "symmetrical"))
  con <- dbConnect(duckdb(), dbdir = "data/results.duckdb")
  dbWriteTable(con, "bf_decision_threshold", results, append = TRUE)
  dbDisconnect(con)
}
