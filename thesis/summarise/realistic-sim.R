# Summaries for the realistic (Cauchy prior) optional stopping simulation of the thesis.
# Queries data/hacking-bayes.duckdb and writes thesis/summaries/; plotted by thesis/plots/realistic-sim.R.
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("shared/summary-io.R")
# Initialise params
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)

db_file <- "data/hacking-bayes.duckdb"

# Summary: average BF per mu and n. Rows are grouped by mu in the order
# returned by SELECT DISTINCT mu, which determines the line colors.
cauchy_summary <- function(summary_path = NULL) {
  # read the simulation results from the file
  con <- dbConnect(duckdb(), db_file)
  # get mus
  mus <- unlist(dbGetQuery(con, "SELECT DISTINCT mu FROM cauchy_prior ORDER BY mu"))
  # average BF over repetitions
  results_df <- dbGetQuery(con, "SELECT mu, n, AVG(bf) AS bf FROM cauchy_prior GROUP BY n, mu ORDER BY mu, n")
  dbDisconnect(con)
  summary_df <- results_df[order(match(results_df$mu, mus)), ]
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "cauchy")

  summary_df
}

# Summary: decision probability for H0 per mu for every (r, bf_crit, trial range)
# combination shown in the overview plot.
realistic_sim_overview_summary <- function(summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)
  settings <- data.frame(
    r = c(r_vals[2], r_vals[2], r_vals[2], r_vals[1], r_vals[3], r_vals[2]),
    bf_crit = c(BF_crits[1], BF_crits[2], BF_crits[3], BF_crits[1], BF_crits[1], BF_crits[1]),
    trial_start = c(2, 2, 2, 2, 2, 20),
    trial_end = c(100000, 100000, 100000, 100000, 100000, 200)
  )
  summary_rows <- vector("list", nrow(settings))
  for (i in seq_len(nrow(settings))) {
    s <- settings[i, ]
    probs <- if (s$trial_start == 2) {
      dbGetQuery(
        con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
        list(s$r, s$bf_crit)
      )
    } else {
      dbGetQuery(
        con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 20 AND trial_end = 200
                GROUP BY mu
                ORDER BY mu",
        list(s$r, s$bf_crit)
      )
    }
    summary_rows[[i]] <- data.frame(
      r = s$r, bf_crit = s$bf_crit,
      trial_start = s$trial_start, trial_end = s$trial_end,
      mu = probs$mu, prob = probs$prob
    )
  }
  dbDisconnect(con)
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "overview")

  summary_df
}

# Plot three BF01 curves for mus and decision probability for H0 given mu
# Summary: stop counts of every repetition per mu and decision
# (0 = H0, 1 = H1, 2 = indecisive) for the given bf_crit and r_val.
realistic_sim_histograms_summary <- function(bf_crit, mus, r_val, summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)
  summary_rows <- list()
  for (mu in mus) {
    h0_data <- dbGetQuery(
      con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND trial_start = 2 AND trial_end = 100000
                    AND ABS(? - mu) < 1e-6
                    AND decision = 0",
      list(r_val, bf_crit, mu)
    )
    h1_data <- dbGetQuery(
      con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND trial_start = 2 AND trial_end = 100000
                    AND ABS(? - mu) < 1e-6
                    AND decision = 1",
      list(r_val, bf_crit, mu)
    )
    indecisive_data <- dbGetQuery(
      con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND trial_start = 2 AND trial_end = 100000
                    AND ABS(? - mu) < 1e-6
                    AND decision = 2",
      list(r_val, bf_crit, mu)
    )
    for (decision in 0:2) {
      stop_counts <- list(h0_data, h1_data, indecisive_data)[[decision + 1]][["stop_count"]]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        bf_crit = rep(bf_crit, length(stop_counts)), r_val = rep(r_val, length(stop_counts)),
        mu = rep(mu, length(stop_counts)), decision = rep(decision, length(stop_counts)),
        stop_count = stop_counts
      )
    }
  }
  dbDisconnect(con)
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "histogram")

  summary_df
}

# Summary: number of simulations and decision probability for H0 per mu.
realistic_sim_decision_prob_curve_summary <- function(bf_crit, r_val, summary_path = NULL) {

  con <- dbConnect(duckdb(), db_file)

  long_range <- dbGetQuery(
    con,
    "
    SELECT
      mu,
      COUNT(*) AS n_total,
      COUNT(CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
    FROM cauchy_sym
    WHERE ABS(? - r) < 1e-6
      AND bf_crit = ?
      AND trial_start = 2
      AND trial_end = 100000
    GROUP BY mu
    ORDER BY mu
    ",
    list(r_val, bf_crit)
  )
  dbDisconnect(con)
  summary_df <- data.frame(bf_crit = bf_crit, r_val = r_val, long_range)

  write_summary(summary_df, summary_path, "decision probability")

  summary_df
}

invisible(realistic_sim_decision_prob_curve_summary(
  bf_crit = 3, r_val = 1 / sqrt(2),
  summary_path = "thesis/summaries/realistic-sim-decision-prob-bf-crit2-3-r-0.707.txt"
))
invisible(cauchy_summary(summary_path = "thesis/summaries/cauchy-prior-simulation.txt"))
invisible(realistic_sim_overview_summary(summary_path = "thesis/summaries/realistic-sym-decision-prob.txt"))
invisible(realistic_sim_histograms_summary(
  bf_crit = 3, c(0.1, 0.5, 0.8), r_val = 1 / sqrt(2),
  summary_path = "thesis/summaries/realistic-sim-bf-crit-3-r-0.707.txt"
))
