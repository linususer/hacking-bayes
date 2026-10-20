# Summaries of the decision probability for H0 by BF_crit and r (thesis, symmetrical rule).
# Queries data/hacking-bayes.duckdb and writes thesis/summaries/; plotted by thesis/plots/decision-prob-bf-crit.R.
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("shared/summary-io.R")

# Initialize variables
big_sim_mus <- c(0.1, 0.2, 0.5, 0.7)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

db_file <- "data/hacking-bayes.duckdb"

# Load results # does not work because file is too big
# sym <- read.csv("data/realistic_sym_simulation.csv")
# Plot the decision probability for H0 given mu
# x-axis: BF_crit, y-axis: P(H0), color: r
# Summary: decision probability for H0 per mu, r and bf_crit.
# Used by plotProbWithBFCritSummary() and plotProbWithBFCritBoundsSummary().
probWithBFCritSummary <- function(mus, summary_path = NULL) {
    con <- dbConnect(duckdb(), db_file, read_only = TRUE)
    summary_rows <- list()
    for (mu in mus) {
        r1_prob <- sapply(BF_crits, function(bf_crit) {
            as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((0.5 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND mu = ?
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit))) / repetitions
        })

        r2_prob <- sapply(BF_crits, function(bf_crit) {
            as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND mu = ?
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit))) / repetitions
        })
        r3_prob <- sapply(BF_crits, function(bf_crit) {
            as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((2 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND mu = ?
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit))) / repetitions
        })
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            mu = mu,
            r = rep(c(0.5, 1, 2) / sqrt(2), each = length(BF_crits)),
            bf_crit = rep(BF_crits, 3),
            prob = c(r1_prob, r2_prob, r3_prob)
        )
    }
    dbDisconnect(con)
    summary_df <- do.call(rbind, summary_rows)
    rownames(summary_df) <- NULL

    write_summary(summary_df, summary_path, "bf_crit")

    summary_df
}

# Summary: decision probability for H0 per mu and bf_crit (r = 1 / sqrt(2)).
decisionProbBFCritsSummary <- function(big_sim_mus, bf_crit, summary_path = NULL) {
    con <- dbConnect(duckdb(), db_file, read_only = TRUE)
    bf1_prob <- sapply(big_sim_mus, function(mu) {
        dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND ROUND(mu, 6) = ROUND(?, 6)
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit[1])) / repetitions
    })
    bf2_prob <- sapply(big_sim_mus, function(mu) {
        dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND ROUND(mu, 6) = ROUND(?, 6)
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit[2])) / repetitions
    })
    bf3_prob <- sapply(big_sim_mus, function(mu) {
        dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND ROUND(mu, 6) = ROUND(?, 6)
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit[3])) / repetitions
    })
    dbDisconnect(con)
    summary_df <- data.frame(
        mu = rep(big_sim_mus, 3),
        bf_crit = rep(bf_crit, each = length(big_sim_mus)),
        prob = c(as.numeric(unlist(bf1_prob)), as.numeric(unlist(bf2_prob)), as.numeric(unlist(bf3_prob)))
    )

    write_summary(summary_df, summary_path, "bf_crits")

    summary_df
}

invisible(probWithBFCritSummary(
    big_sim_mus,
    summary_path = "thesis/summaries/realistic-sym-decision-prob-by-bf_crit.txt"
))
invisible(decisionProbBFCritsSummary(
    seq(0, 1, 0.01), c(3, 6, 10),
    summary_path = "thesis/summaries/realistic-sym-decision-prob-bf-crits.txt"
))
