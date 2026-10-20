# Summaries for the idealised optional stopping simulation (table bf_decision_threshold).
# Queries data/hacking-bayes.duckdb and writes thesis/summaries/; plotted by thesis/plots/optional-stopping.R.
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
library(duckdb)
source("shared/summary-io.R")
#################
### LOAD DATA ###
#################
con <- dbConnect(duckdb(), "data/hacking-bayes.duckdb")
sym <- dbGetQuery(con, "SELECT *
                        FROM bf_decision_threshold
                        WHERE symmetrical = TRUE
                        AND trial_start = 1")
asym <- dbGetQuery(con, "SELECT *
                        FROM bf_decision_threshold
                        WHERE symmetrical = FALSE
                        AND trial_start = 1")
dbDisconnect(con)
##############
### CONFIG ###
##############
big_sim_mus <- seq(0, 1, 0.01)
repetitions <- 20000

##########################
### PLOTTING FUNCTIONS ###
##########################
# Summary: stop counts of every H0 (decision 0) and H1 (decision 1) decision for each mu.
# df_name ("asym" / "sym") is the name of the data passed in and ends up in the file names.
sim_histograms_summary <- function(df, mus, summary_path = NULL) {
  # plot histograms for h0 and h1 decisions
  h0_data <- df[df[, 1] == 0, , ]
  h1_data <- df[df[, 1] == 1, , ]
  df_name <- deparse(substitute(df))
  summary_rows <- list()
  for (mu in mus) {
    h0_counts <- h0_data[h0_data[, 3] == mu, 2]
    h1_counts <- h1_data[h1_data[, 3] == mu, 2]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      df_name = df_name, mu = mu,
      decision = c(rep(0, length(h0_counts)), rep(1, length(h1_counts))),
      stop_count = c(h0_counts, h1_counts)
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "histogram")

  summary_df
}

#################################
### PLOT DECISION PROBABILITY ###
#################################
# Summary: decision probability for H0 per mu for the asymmetrical and/or symmetrical rule.
# case: "asymmetrical", "symmetrical", "both"
decision_prob_curve_summary <- function(case = "both", summary_path = NULL) {
  summary_rows <- list()
  # simulation asymmetrical
  if (case == "asymmetrical" || case == "both") {
    h0_prob <- c()
    h0_data <- asym[asym[, 1] == 0, , ]
    is_ph0_smaller_50 <- FALSE
    for (mu in big_sim_mus) {
      h0_count <- length(h0_data[h0_data[, 3] == mu, 2])
      h0_prob <- c(h0_prob, h0_count / repetitions)
      if (!is_ph0_smaller_50 && h0_count / repetitions < 0.5) {
        is_ph0_smaller_50 <- TRUE
        print(paste("P(H0) > 0.5 for mu =", mu))
      }
    }
    summary_rows[[length(summary_rows) + 1]] <- data.frame(case = case, rule = "asymmetrical", mu = big_sim_mus, h0_prob = h0_prob)
  }
  # simulation symmetrical
  if (case == "symmetrical" || case == "both") {
    h0_prob <- c()
    h0_data <- sym[sym[, 1] == 0, , ]
    is_ph0_smaller_50 <- FALSE
    for (mu in big_sim_mus) {
      h0_count <- length(h0_data[h0_data[, 3] == mu, 2])
      h0_prob <- c(h0_prob, h0_count / repetitions)
      if (!is_ph0_smaller_50 && h0_count / repetitions < 0.5) {
        is_ph0_smaller_50 <- TRUE
        print(paste("P(H0) > 0.5 for mu =", mu))
      }
    }
    summary_rows[[length(summary_rows) + 1]] <- data.frame(case = case, rule = "symmetrical", mu = big_sim_mus, h0_prob = h0_prob)
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  write_summary(summary_df, summary_path, "decision probability")

  summary_df
}

invisible(sim_histograms_summary(asym, c(0.1, 0.5, 0.8), summary_path = "thesis/summaries/asym-both.txt"))
invisible(sim_histograms_summary(sym, c(0.1, 0.2, 0.5, 0.8, 1), summary_path = "thesis/summaries/sym-both.txt"))
for (case in c("asymmetrical", "symmetrical", "both")) {
  decision_prob_curve_summary(
    case = case,
    summary_path = sprintf("thesis/summaries/%s-decision-probability.txt", case)
  )
}
