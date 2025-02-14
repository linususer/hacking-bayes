setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# Load required libraries
library(BayesFactor)
library(data.table)
library(duckdb)
library(foreach)
library(doParallel)

#####################
##### FUNCTIONS #####
#####################

# Calculate Bayes Factor of two hypotheses with binomial distribution density functions
# x: vector of observations
# p1: probability of hypothesis 0 (null hyptoheis)
# p2: probability of hypothesis 1 (alternative hypothesis)
binomBF <- function(x, p0, p1) {
    dbinom(x = sum(x), size = length(x), p = p1) /
        dbinom(x = sum(x), size = length(x), p = p0)
}
probComp <- function(p0, p1, weight) {
    p0 * weight + p1 * (1 - weight)
}

# Function to perform the simulation
# p_true: true probability of binomial distributed data
# p0: probability of hypothesis 0 (null hypothesis)
# p1: probability of hypothesis 1 (alternative hypothesis)
# bf_crit: critical Bayes Factor
# repetitions: number of repetitions
# return a data frame with the results
simulate_bayes_factor <- function(p_true, p0, p1, bf_crit, repetitions, trial_start = 1, trial_end = 100000, sym = FALSE, weight = 0.5) {
    start <- Sys.time()
    results <- vector("list", length = repetitions)
    # Simulation
    for (i in 1:repetitions) {
        x <- NULL
        x <- rbinom(n = 1, size = 1, p_true)
        bf <- binomBF(x, p0, p1)
        stop_count <- 1
        # Asymmetrical Optional Stopping
        if (!sym) {
            while (is.na(bf) || is.null(bf) || (bf > (1 / bf_crit) && stop_count < trial_end)) { # sym: && bf < bf_crit)) {
                x <- c(x, rbinom(n = 1, size = 1, p_true))
                bf <- binomBF(x, p0, p1)
                stop_count <- stop_count + 1
            }
        }
        # Symmetrical Optional Stopping
        else {
            while (is.na(bf) || is.null(bf) || (bf > (1 / bf_crit) && bf < bf_crit && stop_count < trial_end)) {
                x <- c(x, rbinom(n = 1, size = 1, p_true))
                bf <- binomBF(x, p0, p1)
                stop_count <- stop_count + 1
            }
        }
        decision <- if (bf < 1 / bf_crit) {
            0
        } else if (bf > bf_crit) {
            1
        } else {
            2
        }
        results[[i]] <- list(decision, stop_count, bf_crit, p0, p1, p_true, sym, trial_start, trial_end)
    }
    results <- rbindlist(results)
    setnames(results, c("decision", "stop_count", "bf_crit", "p0", "p1", "p_true", "symmetrical", "trial_start", "trial_end"))
    end <- Sys.time()
    print(paste("Simulation done for p_true =", p_true, "with P0 =", p0, "and P1 =", p1, "in", end - start))
    results
}
#####################
### CONFIGURATION ###
#####################

# Replicating Fig 2 from Sanborn (2014)
prob_heads <- 0.75
prob_tails <- 0.25
weight <- 0.5
prob_comp <- probComp(prob_heads, prob_tails, weight)
probs_true <- seq(prob_tails, prob_heads, by = 0.01)
repetitions <- 20000
bf_crit <- 10 # c(3, 6, 10)
db_file <- "data/results.duckdb"

registerDoParallel(cores = detectCores() - 1)
# Create a connection to the DuckDB database
con <- dbConnect(duckdb(), db_file)
dbExecute(con, "DROP TABLE IF EXISTS sanborn_probs_replication")
dbExecute(con, "CREATE TABLE sanborn_probs_replication (decision INTEGER, stop_count INTEGER, bf_crit DOUBLE,
                p0 DOUBLE, p1 DOUBLE, p_true DOUBLE, symmetrical BOOLEAN, trial_start INTEGER, trial_end INTEGER)")

# Simulate Optional Stopping for different true probabilities between [1, 100000]
results <- rbindlist(foreach(prob_true = probs_true) %dopar% {
    simulate_bayes_factor(prob_true, prob_heads, prob_tails, bf_crit, repetitions, sym = TRUE, weight = weight)
})
dbWriteTable(con, "sanborn_probs_replication", results, append = TRUE)

dbDisconnect(con)


# Simulate Optional Stopping for different true probabilities between [20, 200]
con <- dbConnect(duckdb(), db_file)
results <- rbindlist(foreach(prob_true = probs_true) %dopar% {
    simulate_bayes_factor(prob_true, prob_heads, prob_tails, bf_crit, repetitions, sym = TRUE, weight = weight, trial_start = 20, trial_end = 200)
})
dbWriteTable(con, "sanborn_probs_replication", results, append = TRUE)
dbDisconnect(con)
stopImplicitCluster()
