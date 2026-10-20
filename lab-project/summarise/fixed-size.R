# Summaries comparing fixed sample size designs with optional stopping (lab project).
# Queries data/hacking-bayes.duckdb and writes lab-project/summaries/; plotted by lab-project/plots/fixed-size.R.
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("shared/queries.R")
source("shared/summary-io.R")
# Initialise params
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3)

db_file <- "data/hacking-bayes.duckdb"

con <- dbConnect(duckdb(), db_file)
fixed_opt_avg_data <- get_fixed_stop_avg(con, BF_crits[1], r_vals[2])
fixed_weighted_sum_data <- get_fixed_same_dist(con, BF_crits[1], r_vals[2])
fixed_max_data <- get_fixed_max(con, BF_crits[1], r_vals[2])
dbDisconnect(con)

# Summary: for each mu and BF_crit the fixed size decision probabilities per
# trial_count ("fixed_size"), the optional stopping decision probabilities and
# counts per stop_count ("opt_stop") and the STOP AVG, SAME DIST and MAX values
# taken from the given data frames ("stop_avg", "same_dist", "max").
realistic_sim_fixed_size_summary <- function(fixed_opt_avg_data, fixed_weighted_sum_data, fixed_max_data,
                                             summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)
  summary_rows <- list()
  for (mu in c(0.5)) { # c(0.01,0.02,0.06,0.05, 0.1,0.2,0.5)
    for (BF_crit in BF_crits) {
      fixed_size <- dbGetQuery(
        con, "SELECT trial_count, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym_fixed_size
                WHERE ABS(? - mu) < 1e-6 AND ABS(? - r) < 1e-6 AND bf_crit = ?
                GROUP BY trial_count
                ORDER BY trial_count",
        list(mu, r_vals[2], BF_crit)
      )
      opt_stop <- dbGetQuery(
        con, "SELECT stop_count, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob,
                COUNT(stop_count) AS count
                FROM cauchy_sym
                WHERE ABS(? - mu) < 1e-6 AND ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY stop_count
                ORDER BY stop_count",
        list(mu, r_vals[2], BF_crit)
      )
      stop_avg <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu) < 1e-6, ]
      same_dist <- fixed_weighted_sum_data[abs(fixed_weighted_sum_data$delta - mu) < 1e-6, ]
      max_rows <- fixed_max_data[abs(fixed_max_data$mu - mu) < 1e-6, ]
      summary_rows[[length(summary_rows) + 1]] <- rbind(
        data.frame(mu = mu, bf_crit = BF_crit, series = rep("fixed_size", nrow(fixed_size)),
                   n = fixed_size$trial_count, prob = fixed_size$prob, count = NA),
        data.frame(mu = mu, bf_crit = BF_crit, series = rep("opt_stop", nrow(opt_stop)),
                   n = opt_stop$stop_count, prob = opt_stop$prob, count = opt_stop$count),
        data.frame(mu = mu, bf_crit = BF_crit, series = rep("stop_avg", nrow(stop_avg)),
                   n = stop_avg$avg_stop, prob = stop_avg$prob, count = NA),
        data.frame(mu = mu, bf_crit = BF_crit, series = rep("same_dist", nrow(same_dist)),
                   n = NA, prob = same_dist$prob, count = NA),
        data.frame(mu = mu, bf_crit = BF_crit, series = rep("max", nrow(max_rows)),
                   n = max_rows$trial_count, prob = max_rows$prob, count = NA)
      )
    }
  }
  dbDisconnect(con)
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "fixed size")

  summary_df
}

# Summary: for every mu the fixed size maximum of the decision probability for H0
# with its trial_count ("fixed_max") and the optional stopping decision probability
# with its mean stop count ("opt_stop").
realistic_sim_fixed_max_summary <- function(bf_crit, r, summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)
  fixed_max <- dbGetQuery(
    con, "WITH prob_table AS (
                    SELECT MU, trial_count, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                    FROM cauchy_sym_fixed_size
                    WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                    GROUP BY mu, trial_count
                    ORDER BY mu, trial_count
                )
                SELECT *
                  FROM prob_table p
                  WHERE prob = (
                    SELECT MAX(prob)
                    FROM prob_table p2
                    WHERE p2.mu = p.mu
                  )
                    ORDER BY mu, trial_count",
    list(r, bf_crit)
  )
  opt_stop <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob,
                AVG(stop_count) AS mean_count
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r, bf_crit)
  )
  dbDisconnect(con)
  summary_df <- rbind(
    data.frame(bf_crit = bf_crit, r = r, series = rep("fixed_max", nrow(fixed_max)),
               mu = fixed_max[["mu"]], n = fixed_max[["trial_count"]], prob = fixed_max[["prob"]]),
    data.frame(bf_crit = bf_crit, r = r, series = rep("opt_stop", nrow(opt_stop)),
               mu = opt_stop[["mu"]], n = opt_stop[["mean_count"]], prob = opt_stop[["prob"]])
  )

  write_summary(summary_df, summary_path, "fixed max")

  summary_df
}

# Summary: decision probability per effect size for optional stopping ("opt_stop"),
# MAX ("fixed_max"), STOP AVG ("stop_avg") and SAME DIST ("same_dist").
realistic_sim_combined_summary <- function(bf_crit, r_val, db_file = "data/hacking-bayes.duckdb", decision = 0,
                                           summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)

  opt_stop <- get_opt_stop(con, bf_crit, r_val, decision)
  fixed_max <- get_fixed_max(con, bf_crit, r_val, decision)
  fixed_opt_avg <- get_fixed_stop_avg(con, bf_crit, r_val, decision)
  weighted_results <- get_fixed_same_dist(con, bf_crit, r_val, decision)

  dbDisconnect(con)

  summary_df <- rbind(
    data.frame(bf_crit = bf_crit, r_val = r_val, decision = decision, series = rep("opt_stop", nrow(opt_stop)),
               delta = opt_stop[["mu"]], prob = opt_stop[["prob"]]),
    data.frame(bf_crit = bf_crit, r_val = r_val, decision = decision, series = rep("fixed_max", nrow(fixed_max)),
               delta = fixed_max[["mu"]], prob = fixed_max[["prob"]]),
    data.frame(bf_crit = bf_crit, r_val = r_val, decision = decision, series = rep("stop_avg", nrow(fixed_opt_avg)),
               delta = fixed_opt_avg[["delta"]], prob = fixed_opt_avg[["prob"]]),
    data.frame(bf_crit = bf_crit, r_val = r_val, decision = decision, series = rep("same_dist", nrow(weighted_results)),
               delta = weighted_results[["delta"]], prob = weighted_results[["prob"]])
  )

  write_summary(summary_df, summary_path, "combined")

  summary_df
}

# Summary: decision probabilities per effect size for every decision
# (0 = H0, 1 = H1, 2 = indecisive) for MAX ("fixed_max", trial_count chosen for
# decision 0), optional stopping ("opt_stop"), STOP AVG ("stop_avg") and SAME DIST ("same_dist").
decision_sim_fixed_size_summary <- function(bf_crit, r_val, summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)
  fixed_max_d0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)
  summary_rows <- list()
  for (decision in 0:2) {
    if (decision == 0) {
      fixed_max <- fixed_max_d0
    } else {
      fixed_max <- get_prob_by_trial(con, bf_crit, r_val, decision, fixed_max_df = fixed_max_d0)
    }
    opt_stop <- get_opt_stop(con, bf_crit, r_val, decision)
    fixed_opt_avg <- get_fixed_stop_avg(con, bf_crit, r_val, decision)
    weighted_results <- get_fixed_same_dist(con, bf_crit, r_val, decision)
    summary_rows[[length(summary_rows) + 1]] <- rbind(
      data.frame(decision = decision, series = rep("fixed_max", nrow(fixed_max)),
                 delta = fixed_max[["mu"]], trial_count = fixed_max[["trial_count"]], prob = fixed_max[["prob"]]),
      data.frame(decision = decision, series = rep("opt_stop", nrow(opt_stop)),
                 delta = opt_stop[["mu"]], trial_count = rep(NA, nrow(opt_stop)), prob = opt_stop[["prob"]]),
      data.frame(decision = decision, series = rep("stop_avg", nrow(fixed_opt_avg)),
                 delta = fixed_opt_avg[["delta"]], trial_count = rep(NA, nrow(fixed_opt_avg)), prob = fixed_opt_avg[["prob"]]),
      data.frame(decision = decision, series = rep("same_dist", nrow(weighted_results)),
                 delta = weighted_results[["delta"]], trial_count = rep(NA, nrow(weighted_results)), prob = weighted_results[["prob"]])
    )
  }
  dbDisconnect(con)
  summary_df <- data.frame(bf_crit = bf_crit, r_val = r_val, do.call(rbind, summary_rows))
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "all decisions")

  summary_df
}

# Summary: for each mu the fixed size decision probabilities for H0 per trial_count
# ("fixed_size"), the optional stopping decision probabilities and counts per
# stop_count ("opt_stop") and the STOP AVG value from fixed_opt_avg_data ("stop_avg").
realistic_sim_fixed_size_hill_summary <- function(fixed_opt_avg_data, summary_path = NULL) {
  con <- dbConnect(duckdb(), db_file)
  mu_vals <- c(0.01, 0.02, 0.06)
  BF_crit <- 10
  summary_rows <- list()
  for (mu in mu_vals) {
    fixed_size <- dbGetQuery(
      con, "SELECT trial_count, COUNT(CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
              FROM cauchy_sym_fixed_size
              WHERE ABS(? - mu) < 1e-6 AND ABS(? - r) < 1e-6 AND bf_crit = ?
              GROUP BY trial_count
              ORDER BY trial_count",
      list(mu, r_vals[2], BF_crit)
    )
    opt_stop <- dbGetQuery(
      con, "SELECT stop_count, COUNT(CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob,
              COUNT(stop_count) AS count
              FROM cauchy_sym
              WHERE ABS(? - mu) < 1e-6 AND ABS(? - r) < 1e-6 AND bf_crit = ?
              AND trial_start = 2 AND trial_end = 100000
              GROUP BY stop_count
              ORDER BY stop_count",
      list(mu, r_vals[2], BF_crit)
    )
    stop_avg <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu) < 1e-6, ]
    summary_rows[[length(summary_rows) + 1]] <- rbind(
      data.frame(mu = mu, series = rep("fixed_size", nrow(fixed_size)),
                 n = fixed_size$trial_count, prob = fixed_size$prob, count = NA),
      data.frame(mu = mu, series = rep("opt_stop", nrow(opt_stop)),
                 n = opt_stop$stop_count, prob = opt_stop$prob, count = opt_stop$count),
      data.frame(mu = mu, series = rep("stop_avg", nrow(stop_avg)),
                 n = NA, prob = stop_avg$prob, count = NA)
    )
  }
  dbDisconnect(con)
  summary_df <- data.frame(bf_crit = BF_crit, do.call(rbind, summary_rows))
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "fixed size hill")

  summary_df
}

# Summary: mirrored decision probability for H0 per effect size, weighted with
# the Cauchy prior density (one row per plotted point).
prior_distr_summary <- function(con, bf_crit, r_val, df_type = "opt_stop", summary_path = NULL) {
  switch(df_type,
    "opt_stop" = {
      df <- get_opt_stop(con, bf_crit, r_val, decision = 0)
      target_delta <- 0.01
      delta_or_mu <- "mu"
    },
    "fixed_max" = {
      df <- get_fixed_max(con, bf_crit, r_val, decision = 0)
      target_delta <- 0.05
      delta_or_mu <- "mu"
    },
    "fixed_opt_avg" = {
      df <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 0)
      target_delta <- 0.05
      delta_or_mu <- "delta"
    },
    "weighted_sum" = {
      df <- get_fixed_same_dist(con, bf_crit, r_val, decision = 0)
      target_delta <- 0.05
      delta_or_mu <- "delta"
    },
    stop("Invalid df")
  )

  # 1. Remove duplicates and ensure the column name stays consistent
  # We use a dynamic formula so it works for both 'mu' and 'delta'
  formula_str <- as.formula(paste("prob ~", delta_or_mu))
  df <- aggregate(formula_str, data = df, FUN = mean)

  # 2. Filter for equal spaces using target_delta
  df <- df[abs(df[[delta_or_mu]] / target_delta - round(df[[delta_or_mu]] / target_delta)) < 1e-9, ]

  # 3. Calculate weights and densities
  deltas <- df[[delta_or_mu]][df[[delta_or_mu]] > 0]
  delta_0 <- df[[delta_or_mu]][df[[delta_or_mu]] == 0]

  weights <- (target_delta / (pi * r_val * (1 + (deltas / r_val)^2)))
  weights_0 <- (target_delta / (pi * r_val * (1 + (delta_0 / r_val)^2)))

  probs <- df$prob[df[[delta_or_mu]] > 0] * weights
  prob_0 <- df$prob[df[[delta_or_mu]] == 0] * weights_0

  # 4. Mirror the distribution
  deltas_plot <- c(-rev(deltas), delta_0, deltas)
  probs_plot <- c(rev(probs), prob_0, probs)

  print(data.frame(delta = deltas_plot, prob = probs_plot))

  summary_df <- data.frame(
    bf_crit = bf_crit, r_val = r_val, df_type = df_type,
    delta = deltas_plot, prob = probs_plot
  )

  write_summary(summary_df, summary_path, "prior distribution")

  summary_df
}

invisible(realistic_sim_fixed_size_summary(
  fixed_opt_avg_data, fixed_weighted_sum_data, fixed_max_data,
  summary_path = "lab-project/summaries/fixed-decision-STOP_AVG.txt"
))
invisible(realistic_sim_fixed_max_summary(
  3, r_vals[2],
  summary_path = "lab-project/summaries/realistic-sim-fixed-max-bf3.txt"
))
invisible(realistic_sim_combined_summary(
  3, r_vals[3],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-1.414.txt"
))
invisible(realistic_sim_combined_summary(
  3, r_vals[2],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-0.707.txt"
))
invisible(realistic_sim_combined_summary(
  3, r_vals[1],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-0.354.txt"
))
invisible(realistic_sim_combined_summary(
  6, r_vals[2],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-6-r-0.707.txt"
))
invisible(realistic_sim_combined_summary(
  10, r_vals[2],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-10-r-0.707.txt"
))
invisible(decision_sim_fixed_size_summary(
  3, r_vals[2],
  summary_path = "lab-project/summaries/realistic-sim-fixed-size-all-decisions-bf-crit-3-r-0.707.txt"
))
invisible(realistic_sim_fixed_size_hill_summary(
  fixed_opt_avg_data,
  summary_path = "lab-project/summaries/fixed-decision-prob-delta-compare.txt"
))

con <- dbConnect(duckdb(), db_file)
invisible(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "opt_stop",
  summary_path = "lab-project/summaries/prior-distribution_opt_stop.txt"
))
invisible(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "fixed_max",
  summary_path = "lab-project/summaries/prior-distribution_fixed_max.txt"
))
invisible(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "fixed_opt_avg",
  summary_path = "lab-project/summaries/prior-distribution_fixed_opt_avg.txt"
))
invisible(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "weighted_sum",
  summary_path = "lab-project/summaries/prior-distribution_weighted_sum.txt"
))
dbDisconnect(con)
