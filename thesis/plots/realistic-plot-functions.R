# Set working directory
setwd(".")

rm(list = ls())
gc()
library(duckdb)
source("shared/define_colors.R")
# Initialise params
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

db_file <- "data/hacking-bayes.duckdb"
table_name_opt_stop <- "realistic_sym_simulation"
table_name_cauchy_fun <- "cauchy_prior"


#######################
####### PLOTTING ######
#######################

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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved cauchy summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by cauchy_summary() or a path
# to the summary .txt it wrote.
plot_cauchy_summary <- function(summary, filename = "thesis/figures/cauchy-prior-simulation.pdf") {
  results_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  mus <- unique(results_df$mu)
  pdf(filename)
  plot(0, 0,
    type = "n", xlab = "n", ylab = bquote(BF[01]),
    xlim = c(1, 100), log = "y", ylim = c(1 / 10, 5),
    main = bquote(BF[01] * " depending on n and " * mu * " with a Cauchy prior")
  )
  for (i in seq_along(mus)) {
    lines(results_df$n[results_df$mu == mus[i]], 1 / (results_df$bf[results_df$mu == mus[i]]),
      col = my_colors[i], lwd = 2
    )
  }
  abline(h = 1, col = "black", lty = 2)
  legend_labels <- sapply(mus, function(mu) {
    bquote(mu == .(mu))
  })
  legend("topright",
    legend = legend_labels,
    col = my_colors, lty = 1, cex = 0.8
  )
  dev.off()
  invisible(results_df)
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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved overview summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by realistic_sim_overview_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_overview_summary <- function(summary,
                                                filename = "thesis/figures/realistic-sym-decision-prob.pdf") {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  get_probs <- function(r, bf_crit, trial_start = 2) {
    summary_df[abs(summary_df$r - r) < 1e-6 & summary_df$bf_crit == bf_crit &
      summary_df$trial_start == trial_start, c("mu", "prob")]
  }
  pdf(filename, width = 12, height = 8)
  par(mfrow = c(2, 3), mar = c(5, 5, 5, 5))
  # Add one space
  plot.new()
  # base case
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = "Base Case",
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(delta),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
  )
  base <- get_probs(r_vals[2], BF_crits[1])
  lines(base[[1]], base[[2]], col = my_colors[1], lwd = 3)
  abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
  # Add one space
  plot.new()
  # Plot BF_crits
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Varying " * BF[crit]),
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(delta),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
  )
  bf1 <- get_probs(r_vals[2], BF_crits[1])
  bf2 <- get_probs(r_vals[2], BF_crits[2])
  bf3 <- get_probs(r_vals[2], BF_crits[3])
  lines(bf1[[1]], bf1[[2]], col = my_colors[1], lwd = 3)
  lines(bf2[[1]], bf2[[2]], col = my_colors[2], lwd = 3)
  lines(bf3[[1]], bf3[[2]], col = my_colors[3], lwd = 3)
  # Add universal bounds
  abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
  abline(h = 1 / BF_crits[2], col = my_colors[2], lty = 2, lwd = 3)
  abline(h = 1 / BF_crits[3], col = my_colors[3], lty = 2, lwd = 3)

  legend("topright",
    legend = c(bquote(BF[crit] == 3), bquote(BF[crit] == 6), bquote(BF[crit] == 10)),
    col = my_colors, lwd = 3
  )
  # Plot the decision probability for H0 given mu
  # for (bf_crit in BF_crits) {
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = "Varying r-scale",
    ylab = "",
    xlab = bquote(delta),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
  )
  r1 <- get_probs(r_vals[1], BF_crits[1])
  r2 <- get_probs(r_vals[2], BF_crits[1])
  r3 <- get_probs(r_vals[3], BF_crits[1])
  lines(r1[[1]], r1[[2]], col = my_colors[2], lwd = 3)
  lines(r2[[1]], r2[[2]], col = my_colors[1], lwd = 3)
  lines(r3[[1]], r3[[2]], col = my_colors[3], lwd = 3)
  abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
  legend("topright",
    legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
    col = c(my_colors[2], my_colors[1], my_colors[3]), lwd = 3
  )
  # Plot 20-200 vs. 2 to 100000
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = "2 to 100000 vs. 20 to 200",
    ylab = "",
    xlab = bquote(delta),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
  )
  short_range <- get_probs(r_vals[2], BF_crits[1], trial_start = 20)
  long_range <- get_probs(r_vals[2], BF_crits[1])
  lines(long_range[[1]], long_range[[2]], col = my_colors[1], lwd = 3)
  lines(short_range[[1]], short_range[[2]], col = my_colors[2], lwd = 3)
  abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
  legend("topright",
    legend = c("2 to 100000", "20 to 200"),
    col = c(my_colors[1], my_colors[2]), lwd = 3
  )
  dev.off()
  # }
  invisible(summary_df)
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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved histogram summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by realistic_sim_histograms_summary()
# or a path to the summary .txt it wrote. Writes one pdf per mu.
plot_realistic_sim_histograms_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  bf_crit <- summary_df$bf_crit[1]
  r_val <- summary_df$r_val[1]
  mus <- unique(summary_df$mu)
  # Three Histograms for mus with H0 and H1 decisions
  custom_breaks <- seq(0, 4000, by = 1)
  r_text <- round(r_val, 3)
  for (mu in mus) {
    h0_data <- summary_df[summary_df$mu == mu & summary_df$decision == 0, "stop_count", drop = FALSE]
    h1_data <- summary_df[summary_df$mu == mu & summary_df$decision == 1, "stop_count", drop = FALSE]
    pdf(paste("thesis/figures/realistic-sim-bf-crit-", bf_crit, "-r-", r_text, "-mu-", mu, ".pdf", sep = ""))
    par(mfrow = c(2, 1))
    par(mar = c(0, 5, 3, 3))
    # print(h0_data %>% select(Stop_Count) %>% pull())
    hist(h0_data[["stop_count"]],
      main = bquote("Simulation with 20000 repetitions for " * mu * " = " * .(mu) * ", " * BF[crit] * " = " * .(bf_crit) * ", " * r * " = " * .(r_text)),
      xlim = c(0, 60), ylim = c(0, 4000),
      xlab = "",
      ylab = bquote("Decision Count (" * H[0] * ")"),
      xaxt = "n", las = 1,
      col = "#59B3E6", breaks = custom_breaks
    )
    legend("topright", legend = c(bquote(H[0]), bquote(H[1])), fill = c("#59B3E6", "#CD1076"))
    par(mar = c(5, 5, 0, 3))
    hist(h1_data[["stop_count"]],
      main = "",
      xlim = c(0, 60), ylim = c(4000, 0),
      xlab = "Stop Count",
      ylab = bquote("Decision Count " * "(" * H[1] * ")"),
      las = 1,
      col = "#CD1076", breaks = custom_breaks
    )
    dev.off()
  }
  invisible(summary_df)
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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved decision probability summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by
# realistic_sim_decision_prob_curve_summary() or a path to the summary .txt it wrote.
plot_realistic_sim_decision_prob_curve_summary <- function(summary) {
  long_range <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  bf_crit <- long_range$bf_crit[1]
  r_val <- long_range$r_val[1]

  r_text <- round(r_val, 3)

  pdf(
    paste0(
      "thesis/figures/realistic-sim-decision-prob-bf-crit2-",
      bf_crit, "-r-", r_text, ".pdf"
    ),
    width = 14, height = 8
  )

  par(mar = c(5, 7, 5, 7))  # extra space for both y-axes

  # total simulations per mu (assumed constant)
  N <- unique(long_range$n_total)

  ## ---- Base plot (probability scale) ----
  plot(
    long_range$mu, long_range$prob,
    type = "l", lwd = 5, col = "black",
    yaxt = "n",
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = bquote("Korrelationsstärke " * delta),
    ylab = "",
    cex.axis = 1.8,
    cex.lab = 2,
    cex.main = 2
  )

  ## ---- Right axis: probability ----
  axis(
    side = 4,
    at = seq(0, 1, by = 0.1),
    labels = paste0(seq(0, 100, by = 10), "%"),
    cex.axis = 1.8
  )
  mtext(
    bquote("Entscheidungswahrscheinlichkeit für Hypothese " * H[0]),
    side = 4, line = 4, cex = 2
  )
  ## ---- Left axis: absolute counts ----
  axis(
    side = 2,
    at = seq(0, 1, by = 0.1),
    labels = seq(0, N, by = N / 10),
    cex.axis = 1.8
  )
  mtext(
    bquote("Entscheidungen für Hypothese " * H[0]),
    side = 2, line = 4, cex = 2
  )

  box()

  dev.off()
  invisible(long_range)
}
decision_prob_summary <- realistic_sim_decision_prob_curve_summary(
  bf_crit = 3, r_val = 1 / sqrt(2),
  summary_path = "thesis/summaries/realistic-sim-decision-prob-bf-crit2-3-r-0.707.txt"
)
plot_realistic_sim_decision_prob_curve_summary(decision_prob_summary)
plot_cauchy_summary(cauchy_summary(summary_path = "thesis/summaries/cauchy-prior-simulation.txt"))
plot_realistic_sim_overview_summary(realistic_sim_overview_summary(
  summary_path = "thesis/summaries/realistic-sym-decision-prob.txt"
))
plot_realistic_sim_histograms_summary(realistic_sim_histograms_summary(
  bf_crit = 3, c(0.1, 0.5, 0.8), r_val = 1 / sqrt(2),
  summary_path = "thesis/summaries/realistic-sim-bf-crit-3-r-0.707.txt"
))
