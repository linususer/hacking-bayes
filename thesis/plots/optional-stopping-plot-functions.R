setwd(".")
# clear workspace
rm(list = ls())
library(duckdb)
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
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
big_sim_mus <- seq(0, 1, 0.01)
fixed_sigma <- 1
sigmas <- c(0.01, 0.1, 1, 5, 10)
max_stop_count <- 500
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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved histogram summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by sim_histograms_summary()
# or a path to the summary .txt it wrote. Writes one pdf per mu.
plot_sim_histograms_summary <- function(summary, xlim_max = 60, ylimh0 = 2500, ylimh1 = 200) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  df_name <- summary_df$df_name[1]
  mus <- unique(summary_df$mu)
  custom_breaks <- seq(0, 4000, by = 1)
  # plot both decisions
  for (mu in mus) {
    pdf(paste("thesis/figures/", df_name, "-both-", mu, ".pdf", sep = ""))
    par(mfrow = c(2, 1))
    par(mar = c(0, 5, 3, 3))
    hist(summary_df$stop_count[summary_df$decision == 0 & summary_df$mu == mu],
      main = bquote("Simulation with 20000 repetitions for " * mu == .(mu)),
      xlim = c(0, xlim_max), ylim = c(0, ylimh0),
      xlab = "",
      ylab = bquote("Decision Count (" * H[0] * ")"),
      xaxt = "n", las = 1,
      col = "#59B3E6", breaks = custom_breaks
    )
    legend("topright", legend = c(bquote(H[0]), bquote(H[1])), fill = c("#59B3E6", "#CD1076"))
    par(mar = c(5, 5, 0, 3))
    hist(summary_df$stop_count[summary_df$decision == 1 & summary_df$mu == mu],
      main = "",
      xlim = c(0, xlim_max), ylim = c(ylimh1, 0),
      xlab = "Stop Count",
      ylab = bquote("Decision Count (" * H[1] * ")"),
      las = 1,
      col = "#CD1076", breaks = custom_breaks
    )
    dev.off()
  }
  invisible(summary_df)
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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved decision probability summary to:", summary_path, "\n")
  }

  summary_df
}

# plot the decision probability for H0 given mu and n_start
# `summary` is either the data.frame returned by decision_prob_curve_summary()
# or a path to the summary .txt it wrote.
# n_start: number of starts for the random walk
plot_decision_prob_curve_summary <- function(summary, n_start) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  case <- summary_df$case[1]
  # plot the decision probability for H0 given mu
  pdf(paste("thesis/figures/", case, if (n_start > 1) {
    paste("-n_start-", n_start, sep = "")
  }, "-decision-probability.pdf", sep = ""))
  main_text <- if (case == "both") {
    " rules"
  } else {
    " rule"
  }
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability for " * .(case) * .(main_text)),
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
  )
  # simulation asymmetrical
  if (case == "asymmetrical" || case == "both") {
    rows <- summary_df[summary_df$rule == "asymmetrical", ]
    lines(rows$mu, rows$h0_prob, col = "orange", lwd = 2)
  }
  # simulation symmetrical
  if (case == "symmetrical" || case == "both") {
    rows <- summary_df[summary_df$rule == "symmetrical", ]
    lines(rows$mu, rows$h0_prob, col = "#FF99FF", lwd = 2)
  }

  # legend
  if (case == "both") {
    legend("topright",
      legend = c("Asymmetrical Simulation", "Symmetrical Simulation"),
      fill = c("orange", "#FF99FF")
    )
  } else if (case == "asymmetrical") {
    legend("topright",
      legend = c("Asymmetrical Simulation"),
      fill = c("orange")
    )
  } else {
    legend("topright",
      legend = c("Symmetrical Simulation"),
      fill = c("#FF99FF")
    )
  }
  dev.off()
  invisible(summary_df)
}

mus <- c(0.1, 0.5, 0.8)
asym_histograms_summary <- sim_histograms_summary(asym, mus, summary_path = "thesis/summaries/asym-both.txt")
plot_sim_histograms_summary(asym_histograms_summary, xlim_max = 60, ylimh0 = 2500, ylimh1 = 2500)
# the committed sym figures use all five config mus and the default H1 axis limit
sym_histograms_summary <- sim_histograms_summary(sym, c(0.1, 0.2, 0.5, 0.8, 1), summary_path = "thesis/summaries/sym-both.txt")
plot_sim_histograms_summary(sym_histograms_summary, xlim_max = 60, ylimh0 = 2500)
for (case in c("asymmetrical", "symmetrical", "both")) {
  plot_decision_prob_curve_summary(decision_prob_curve_summary(
    case = case,
    summary_path = sprintf("thesis/summaries/%s-decision-probability.txt", case)
  ), n_start = 1)
}
