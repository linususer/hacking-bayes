setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# Load required libraries
library(BayesFactor)
library(data.table)
library(duckdb)

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

alternativeProbComp <- function(p0, p1, weight) {
    heads <- rbinom(n = 1, size = 1, weight)
    if (heads == 1) {
        rbinom(n = 1, size = 1, p0)
    } else {
        rbinom(n = 1, size = 1, p1)
    }
}

# Function to perform the simulation
# p_true: true probability of binomial distributed data
# p0: probability of hypothesis 0 (null hypothesis)
# p1: probability of hypothesis 1 (alternative hypothesis)
# bf_crit: critical Bayes Factor
# repetitions: number of repetitions
# return a data frame with the results
simulate_bayes_factor <- function(p_true, p0, p1, bf_crit, repetitions, sym = FALSE, alt_comp = FALSE, weight = 0.5, alt_id = "") {
    start <- Sys.time()
    results <- vector("list", length = repetitions)
    # Simulation
    for (i in 1:repetitions) {
        x <- NULL
        if (alt_comp) {
            x <- alternativeProbComp(p0, p1, weight)
        } else {
            x <- rbinom(n = 1, size = 1, p_true)
        }
        bf <- binomBF(x, p0, p1)
        stop_count <- 1
        # Asymmetrical Optional Stopping
        if (!sym) {
            while (is.na(bf) || is.null(bf) || (bf > (1 / bf_crit) && stop_count < 250)) { # sym: && bf < bf_crit)) {
                if (alt_comp) {
                    x <- c(x, alternativeProbComp(p0, p1, weight))
                } else {
                    x <- c(x, rbinom(n = 1, size = 1, p_true))
                }
                bf <- binomBF(x, p0, p1)
                stop_count <- stop_count + 1
            }
        }
        # Symmetrical Optional Stopping
        else {
            while (is.na(bf) || is.null(bf) || (bf > (1 / bf_crit) && bf < bf_crit)) {
                if (alt_comp) {
                    x <- c(x, alternativeProbComp(p0, p1, weight))
                } else {
                    x <- c(x, rbinom(n = 1, size = 1, p_true))
                }
                bf <- binomBF(x, p0, p1)
                stop_count <- stop_count + 1
            }
        }
        decision <- ifelse(bf < bf_crit, 0, 1)
        results[[i]] <- list(decision, stop_count, bf_crit, p0, p1, p_true, sym, alt_comp, alt_id)
    }
    results <- rbindlist(results)
    setnames(results, c("decision", "stop_count", "bf_crit", "p0", "p1", "p_true", "symmetrical", "alt_comp", "alt_id"))
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
repetitions <- 20000
bf_crits <- 10 # c(3, 6, 10)
db_file <- "data/results.duckdb"
#####################
#### SIMULATION #####
#####################

con <- dbConnect(duckdb(), dbdir = db_file, read_only = FALSE)
dbExecute(con, "DROP TABLE IF EXISTS sanborn_replication")
dbExecute(con, "CREATE TABLE sanborn_replication (decision INTEGER, stop_count INTEGER, bf_crit DOUBLE,
                p0 DOUBLE, p1 DOUBLE, p_true DOUBLE, symmetrical BOOLEAN, alt_comp BOOLEAN, alt_id STRING)")
dbDisconnect(con)

####### ASYM, p_truth = p0 = 0.75, p1 = 0.25########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_heads, p0 = prob_heads, p1 = prob_tails,
        bf_crit = bf_crit, repetitions = repetitions
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### ASYM, p_truth =0.75, p0 = 0.25, p1 = 0.75########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_heads, p0 = prob_tails, p1 = prob_heads,
        bf_crit = bf_crit, repetitions = repetitions
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### ASYM, p_truth = 0.5, p0 = 0.75, p1 = 0.25 ########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_comp, p0 = prob_heads, p1 = prob_tails,
        bf_crit = bf_crit, repetitions = repetitions
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})
####### ASYM, p_truth = 0.5, p0 = 0.25, p1 = 0.75 ########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_comp, p0 = prob_tails, p1 = prob_heads,
        bf_crit = bf_crit, repetitions = repetitions
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### ASYM, p_truth = ALT_COMP, p0 = 0.75, p1 = 0.25 ########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_comp, p0 = prob_heads, p1 = prob_tails,
        bf_crit = bf_crit, repetitions = repetitions, alt_comp = TRUE, alt_id = "asym_comp1"
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### ASYM, p_truth = ALT_COMP, p0 = 0.25, p1 = 0.75 ########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_comp, p0 = prob_tails, p1 = prob_heads,
        bf_crit = bf_crit, repetitions = repetitions, alt_comp = TRUE, alt_id = "asym_comp2"
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### SYM, p_truth = p0 = 0.75, p1 = 0.25 ########
# Perform parallel simulations
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_heads, p0 = prob_heads, p1 = prob_tails,
        bf_crit = bf_crit, repetitions = repetitions, sym = TRUE
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### SYM, p_truth = 0.5, p0 = 0.75, p1 = 0.25 ########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_comp, p0 = prob_heads, p1 = prob_tails,
        bf_crit = bf_crit, repetitions = repetitions, sym = TRUE
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})

####### SYM, p_truth = ALT_COMP, p0 = 0.75, p1 = 0.25 ########
# Perform parallel simulations
sapply(bf_crits, function(bf_crit) {
    result <- simulate_bayes_factor(
        p_true = prob_comp, p0 = prob_heads, p1 = prob_tails,
        bf_crit = bf_crit, repetitions = 20000, sym = TRUE, alt_comp = TRUE, alt_id = "sym_comp1"
    )
    con <- dbConnect(duckdb(), dbdir = db_file)
    dbWriteTable(con, "sanborn_replication", result, append = TRUE)
    dbDisconnect(con)
})
