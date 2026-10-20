# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("shared/define_colors.R")
source("shared/queries.R")
# Initialise params
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

db_file <- "data/hacking-bayes.duckdb"

###############
#### PLOTS ####
################

fixed_big_mus <- c(seq(0, 0.09, 0.01), seq(0.1, 1, 0.05))
BF_crits <- c(3)

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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved fixed size summary to:", summary_path, "\n")
  }

  summary_df
}

# Plot decision probability for 'H0' on y axis and sample size on the x axis
# `summary` is either the data.frame returned by realistic_sim_fixed_size_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_fixed_size_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  # increase the size, including font size
  for (mu in unique(summary_df$mu)) { # c(0.01,0.02,0.06,0.05, 0.1,0.2,0.5)
    pdf(paste("lab-project/figures/report/fixed-decision-STOP_AVG.pdf", sep = ""))
    par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
    plot(0, 0,
      xlim = c(0, 100), ylim = c(0, 1), type = "n",
      main = bquote("STOP AVG, " * delta == .(mu)),
      ylab = bquote("Decision Probability for " * H[0]), xlab = "fixed sample size / stop count n",
      cex.main = 2, cex.lab = 2, cex.axis = 2
    )
    for (BF_crit in unique(summary_df$bf_crit[summary_df$mu == mu])) {
      rows <- summary_df[summary_df$mu == mu & summary_df$bf_crit == BF_crit, ]
      fixed_size <- data.frame(
        trial_count = rows$n[rows$series == "fixed_size"], prob = rows$prob[rows$series == "fixed_size"]
      )
      opt_stop <- data.frame(
        stop_count = rows$n[rows$series == "opt_stop"], prob = rows$prob[rows$series == "opt_stop"],
        count = rows$count[rows$series == "opt_stop"]
      )
      opt_stop_mean <- sum((opt_stop[[1]] * opt_stop[[3]])) / 20000
      opt_stop_sd <- sd(opt_stop[[1]])
      opt_prob_mean <- sum((opt_stop[[2]] * opt_stop[[3]])) / 20000
      opt_prob_sd <- sd(opt_stop[[2]])

      stop_counts_expanded <- rep(opt_stop$stop_count, opt_stop$count)
      max_stop_count <- max(opt_stop$stop_count)

      # Create histogram data with breaks
      hist_data <- hist(stop_counts_expanded, breaks = c(seq(0, max_stop_count, 1)), plot = FALSE)

      # Norm histogram max to 1
      hist_heights <- hist_data$counts / max(hist_data$counts)

      # plot histogram
      for (i in 1:length(hist_data$counts)) {
        rect(
          xleft = hist_data$breaks[i],
          xright = hist_data$breaks[i + 1],
          ybottom = 0,
          ytop = hist_heights[i],
          col = my_colors[7], # adjustcolor(my_colors[7], alpha.f = 0.7),
          border = NA
        )
      }
      # abline(h = opt_prob_mean, col = my_colors[2], lwd = 5, lty = 2)
      lines(fixed_size[[1]], fixed_size[[2]], col = my_colors[1], lwd = 5)
      # abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
      interp_val <- rows$prob[rows$series == "stop_avg"]
      trial_at_interp_val <- rows$n[rows$series == "stop_avg"]
      weighted_val <- rows$prob[rows$series == "same_dist"]
      max_val <- rows$prob[rows$series == "max"]
      trial_at_max_val <- rows$n[rows$series == "max"]

      # abline(h = max_val, col = my_colors[3], lwd = 5, lty = 3)
      # points(trial_at_max_val, max_val, col = my_colors[3], pch = 4, lwd = 4, cex = 2)

      if (length(interp_val) == 1) {
        abline(h = interp_val, col = my_colors[4], lwd = 5, lty = 4)
        points(trial_at_interp_val, interp_val, col = my_colors[4], pch = 4, lwd = 4, cex = 2)
        lines(c(trial_at_interp_val, trial_at_interp_val), c(0, interp_val), col = my_colors[4], lty = 2, lwd = 5)
      }
      # if (length(weighted_val) == 1) {
      #   abline(h = weighted_val, col = my_colors[5], lwd = 5, lty = 6)
      # }
      legend("topright",
        legend = c(
          "fixed n", "stop count distribution",
          # "Optional Stopping",
          # "MAX",
          "STOP AVG"
          # "SAME DIST"
        ),
        col = c(
          my_colors[1], my_colors[7],
          # my_colors[2],
          # my_colors[3],
          my_colors[4]
          # my_colors[5]
        ),
        lwd = 5, lty = c(
          1, 1,
          #  2,
          #  3
          4
          # 6
        )
      )
    }
    dev.off()
  }
  invisible(summary_df)
}
fixed_size_summary <- realistic_sim_fixed_size_summary(
  fixed_opt_avg_data, fixed_weighted_sum_data, fixed_max_data,
  summary_path = "lab-project/summaries/fixed-decision-STOP_AVG.txt"
)
plot_realistic_sim_fixed_size_summary(fixed_size_summary)

# Plot maximum probability of 'H0' for each mu

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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved fixed max summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by realistic_sim_fixed_max_summary()
# or a path to the summary .txt it wrote. Draws into the current device.
plot_realistic_sim_fixed_max_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  fixed_max <- data.frame(
    mu = summary_df$mu[summary_df$series == "fixed_max"],
    trial_count = summary_df$n[summary_df$series == "fixed_max"],
    prob = summary_df$prob[summary_df$series == "fixed_max"]
  )
  opt_stop <- data.frame(
    mu = summary_df$mu[summary_df$series == "opt_stop"],
    prob = summary_df$prob[summary_df$series == "opt_stop"],
    mean_count = summary_df$n[summary_df$series == "opt_stop"]
  )
  # Decision probability for H0 given mu
  # pdf(paste("lab-project/figures/realistic-sim-fixed-max-bf3.pdf", sep = ""))
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability for 'H0'"),
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote("Effect size " * delta)
  )
  # get the first index where the probability is smaller than 0.5
  lines(fixed_max[["mu"]], fixed_max[["prob"]], col = "black", lwd = 2)
  lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd = 2)
  legend("topright", legend = c("fixed probability maximum", "optional stopping probability"), fill = c("black", my_colors[2]))
  # dev.off()
  max_arg <- max_arg <- merge(
    fixed_max[, c("mu", "trial_count")],
    opt_stop[, c("mu", "mean_count")],
    by = "mu",
    suffixes = c("_fixed", "_opt_stop")
  )
  print(max_arg)
  invisible(summary_df)
}
pdf("lab-project/figures/realistic-sim-fixed-max-bf3.pdf")
plot_realistic_sim_fixed_max_summary(realistic_sim_fixed_max_summary(
  3, r_vals[2],
  summary_path = "lab-project/summaries/realistic-sim-fixed-max-bf3.txt"
))
dev.off()


# P('Diff') = P(H0 | n_opt_stop_avg) and P(H0 | n_fixed) for n_fixed ~= n_opt_stop_avg
# get stop_counts from optional stopping results
#  P('Diff')


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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved combined summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by realistic_sim_combined_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_combined_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  bf_crit <- summary_df$bf_crit[1]
  r_val <- summary_df$r_val[1]
  series <- function(name) summary_df[summary_df$series == name, ]
  opt_stop <- data.frame(mu = series("opt_stop")$delta, prob = series("opt_stop")$prob)
  fixed_max <- data.frame(mu = series("fixed_max")$delta, prob = series("fixed_max")$prob)
  fixed_opt_avg <- data.frame(delta = series("stop_avg")$delta, prob = series("stop_avg")$prob)
  weighted_results <- data.frame(delta = series("same_dist")$delta, prob = series("same_dist")$prob)

  # Plot
  pdf(paste("lab-project/figures/report/realistic-fixed-sim-overview-for-bfcrit-", bf_crit, "-r-", round(r_val, 3), ".pdf", sep = ""))
  par(mar = c(5, 6, 5, 5))
  # r_val text
  if (abs(r_val - r_vals[1]) < 1e-6) {
    r_text <- bquote(r == 0.5 / sqrt(2))
  } else if (abs(r_val - r_vals[2]) < 1e-6) {
    r_text <- bquote(r == 1 / sqrt(2))
  } else if (abs(r_val - r_vals[3]) < 1e-6) {
    r_text <- bquote(r == 2 / sqrt(2))
  }
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote(BF[crit] == .(bf_crit) * ", " * .(r_text)),
    ylab = bquote("Decision Probability for " * H[0]),
    xlab = bquote("Effect size " * delta),
    cex.main = 2, cex.lab = 2, cex.axis = 2
  )

  # print values for stop_avg

  print(data.frame(delta = fixed_opt_avg[["delta"]], prob = fixed_opt_avg[["prob"]]))

  lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd = 5, lty = 1)
  lines(fixed_max[["mu"]], fixed_max[["prob"]], col = my_colors[3], lwd = 5, lty = 3)
  lines(fixed_opt_avg[["delta"]], fixed_opt_avg[["prob"]], col = my_colors[4], lwd = 5, lty = 4)
  lines(weighted_results[["delta"]], weighted_results[["prob"]], col = my_colors[5], lwd = 5, lty = 6)

  ### FUTURE ME: THIS IS CLAUDE CODE, DONT HATE ME!
  # Find intersection between opt_stop and fixed_max
  # Interpolate fixed_max probabilities at opt_stop's mu values
  fixed_max_interp <- approx(fixed_max[["mu"]], fixed_max[["prob"]], xout = opt_stop[["mu"]])$y

  # Find where the sign of the difference changes (crossing point)
  diff <- opt_stop[["prob"]] - fixed_max_interp
  sign_change <- which(diff[-1] * diff[-length(diff)] < 0)[1]

  # Linear interpolation to get precise x intersection
  x1 <- opt_stop[["mu"]][sign_change]
  x2 <- opt_stop[["mu"]][sign_change + 1]
  y1 <- diff[sign_change]
  y2 <- diff[sign_change + 1]
  x_intersect <- x1 - y1 * (x2 - x1) / (y2 - y1)
  y_intersect <- approx(fixed_max[["mu"]], fixed_max[["prob"]], xout = x_intersect)$y

  # Draw arrow pointing to intersection (arrow comes from upper-right offset)
  arrows(
    x0 = x_intersect + 0.1, y0 = y_intersect + 0.1,
    x1 = x_intersect + 0.01, y1 = y_intersect + 0.01,
    length = 0.15, lwd = 3, col = "black"
  )
  ### CLAUDE CODE END
  legend("topright",
    legend = c("optional stopping", "MAX", "STOP AVG", "SAME DIST"),
    col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
    lwd = 5, lty = c(1, 3, 4, 6)
  )
  dev.off()
  invisible(summary_df)
}

plot_realistic_sim_combined_summary(realistic_sim_combined_summary(
  3, r_vals[3],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-1.414.txt"
))
plot_realistic_sim_combined_summary(realistic_sim_combined_summary(
  3, r_vals[2],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-0.707.txt"
))
plot_realistic_sim_combined_summary(realistic_sim_combined_summary(
  3, r_vals[1],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-0.354.txt"
))
plot_realistic_sim_combined_summary(realistic_sim_combined_summary(
  6, r_vals[2],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-6-r-0.707.txt"
))
plot_realistic_sim_combined_summary(realistic_sim_combined_summary(
  10, r_vals[2],
  summary_path = "lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-10-r-0.707.txt"
))

# stop_count_to_effect_size <- function(con, bf_crit, r_val) {
#   # Query with rounding and correct GROUP BY
#   stop_count_data <- dbGetQuery(
#     con, "
#       SELECT mu, AVG(stop_count) AS avg_stop_count
#       FROM cauchy_sym
#       WHERE bf_crit = ? AND ABS(r - ?) < 1e-6
#       AND trial_start = 2 AND trial_end = 100000
#       GROUP BY mu
#       ORDER BY mu
#     ",
#     list(bf_crit, r_val)
#   )
#   stop_count_data
# }

# # plot for all rs and for all BF_crits
# plot_stop_count_to_effect_size_bf <- function(bf_crits, r_val, db_file = "data/hacking-bayes.duckdb") {
#   con <- dbConnect(duckdb(), db_file)
#   pdf(paste("lab-project/figures/stop-count-to-effect-size-r-", round(r_val, 3), ".pdf", sep = ""), width = 12, height = 8)
#   plot(0, 0,
#     xlim = c(0, 1), ylim = c(0, 400), type = "n",
#     main = bquote("Stop Count to Effect Size for " * r * " = " * .(r_val)),
#     ylab = bquote("Stop Count"),
#     xlab = bquote("Effect Size " * delta)
#   )
#   for (i in seq_along(bf_crits)) {
#     stop_count_data <- stop_count_to_effect_size(con, bf_crits[i], r_val)
#     # Plot the stop count data
#     lines(stop_count_data$mu, stop_count_data$avg_stop_count, col = my_colors[i], lwd = 2)
#   }
#   dbDisconnect(con)
#   # Add a legend entry for each BF_crit
#   legend("topright",
#     legend = c(bquote(BF[crit] == .(bf_crits[1])), bquote(BF[crit] == .(bf_crits[2])), bquote(BF[crit] == .(bf_crits[3]))),
#     col = c(my_colors[1], my_colors[2], my_colors[3]), lwd = 2, cex = 0.8
#   )
#   dev.off()
# }

# plot for all rs
# plot_stop_count_to_effect_size_r <- function(bf_crit, r_vals, db_file = "data/hacking-bayes.duckdb") {
#   con <- dbConnect(duckdb(), db_file)
#   pdf(paste("lab-project/figures/stop-count-to-effect-size-bf-", bf_crit, ".pdf", sep = ""), width = 12, height = 8)
#   plot(0, 0,
#     xlim = c(0, 1), ylim = c(0, 50), type = "n",
#     main = bquote("Stop Count to Effect Size for " * BF[crit] * " = " * .(bf_crit)),
#     ylab = bquote("Stop Count"),
#     xlab = bquote("Effect Size " * delta)
#   )
#   for (i in seq_along(r_vals)) {
#     stop_count_data <- stop_count_to_effect_size(con, bf_crit, r_vals[i])
#     # Plot the stop count data
#     lines(stop_count_data$mu, stop_count_data$avg_stop_count, col = my_colors[i], lwd = 2)
#   }
#   dbDisconnect(con)
#   # Add a legend entry for each r_val
#   legend("topright",
#     legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
#     col = c(my_colors[1], my_colors[2], my_colors[3]), lwd = 2, cex = 0.8
#   )
#   dev.off()
# }

# plot_stop_count_to_effect_size_bf(c(3, 6, 10), 1 / sqrt(2))
# plot_stop_count_to_effect_size_r(3, c(0.5 / sqrt(2), 1 / sqrt(2), 2 / sqrt(2)))

# plot difference between optional stopping and fixed size for the same mu
# plot_fixed_vs_optional_stopping <- function(bf_crit, r_val, db_file = "data/hacking-bayes.duckdb") {
#   con <- dbConnect(duckdb(), db_file)
#   pdf(paste("lab-project/figures/fixed-vs-optional-stopping-bf-", bf_crit, "-r-", round(r_val, 3), ".pdf", sep = ""), width = 12, height = 8)
#   plot(0, 0,
#     xlim = c(0, 1), ylim = c(0, 1), type = "n",
#     main = bquote("Fixed Size vs Optional Stopping for " * BF[crit] * " = " * .(bf_crit) * ", r = " * .(r_val)),
#     ylab = bquote("Decision Probability for " * H[0]),
#     xlab = bquote("Effect Size " * delta)
#   )

#   fixed_max <- get_fixed_max(con, bf_crit, 2 * r_val)
#   fixed_opt_avg <- get_fixed_stop_avg(con, bf_crit, 2 * r_val)
#   weighted_sum <- get_fixed_same_dist(con, bf_crit, 2 * r_val)
#   opt_stop <- get_opt_stop(con, bf_crit, r_val)

#   lines(fixed_max$mu, fixed_max$prob, col = my_colors[3], lwd = 2, lty = 2)
#   lines(opt_stop$mu, opt_stop$prob, col = my_colors[2], lwd = 2)
#   lines(fixed_opt_avg$delta, fixed_opt_avg$prob, col = my_colors[4], lwd = 2, lty = 4)
#   lines(weighted_sum$delta, weighted_sum$weighted_sum, col = my_colors[5], lwd = 2, lty = 5)

#   legend("topright",
#     legend = c("optional stopping prob", bquote(max[fixed] * " " * n), bquote(bar(n)[os] * "as fixed size"), "same distribution as fixed"),
#     col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
#     lwd = 2, lty = c(1, 2, 4, 5)
#   )

#   dev.off()

#   dbDisconnect(con)
# }
# plot_fixed_vs_optional_stopping(3, r_vals[2])

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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved all decisions summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by decision_sim_fixed_size_summary()
# or a path to the summary .txt it wrote. Writes one pdf per decision.
plot_decision_sim_fixed_size_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  bf_crit <- summary_df$bf_crit[1]
  r_val <- summary_df$r_val[1]
  for (decision in 0:2) {
    pdf_filename <- sprintf("lab-project/figures/report/realistic-sim-fixed-size-all-decisions-bf-crit-%s-r-%.3f-%s-wo-legend.pdf", bf_crit, r_val, decision)
    pdf(pdf_filename)

    decision_text <- switch(as.character(decision),
      "0" = bquote(H[0]),
      "1" = bquote(H[1]),
      "2" = bquote("Indecisive")
    )

    series <- function(name) summary_df[summary_df$decision == decision & summary_df$series == name, ]
    fixed_max <- data.frame(mu = series("fixed_max")$delta, trial_count = series("fixed_max")$trial_count,
                            prob = series("fixed_max")$prob)
    print(fixed_max)
    opt_stop <- data.frame(mu = series("opt_stop")$delta, prob = series("opt_stop")$prob)
    fixed_opt_avg <- data.frame(delta = series("stop_avg")$delta, prob = series("stop_avg")$prob)
    weighted_results <- data.frame(delta = series("same_dist")$delta, prob = series("same_dist")$prob)
    par(mar = c(5, 6, 5, 5))
    plot(0, 0,
      xlim = c(0, 1), ylim = c(0, 1), type = "n",
      main = decision_text,
      ylab = bquote("Decision Probability for " * .(decision_text)),
      xlab = bquote("Effect size " * delta),
      cex.main = 2, cex.lab = 2, cex.axis = 2
    )

    if (nrow(opt_stop) > 0) lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd = 5, lty = 2)
    if (nrow(fixed_max) > 0) lines(fixed_max[["mu"]], fixed_max[["prob"]], col = my_colors[3], lwd = 5, lty = 3)
    if (nrow(fixed_opt_avg) > 0) lines(fixed_opt_avg[["delta"]], fixed_opt_avg[["prob"]], col = my_colors[4], lwd = 5, lty = 4)
    if (nrow(weighted_results) > 0) lines(weighted_results[["delta"]], weighted_results[["prob"]], col = my_colors[5], lwd = 5, lty = 6)

    # legend("topright",
    #   legend = c("optional stopping", "MAX", "STOP AVG", "SAME DIST"),
    #   col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
    #   lwd = 5, lty = c(2, 3, 4, 6)
    # )
    dev.off()
  }
  invisible(summary_df)
}

plot_decision_sim_fixed_size_summary(decision_sim_fixed_size_summary(
  3, r_vals[2],
  summary_path = "lab-project/summaries/realistic-sim-fixed-size-all-decisions-bf-crit-3-r-0.707.txt"
))

BF_crits <- c(10)
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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved fixed size hill summary to:", summary_path, "\n")
  }

  summary_df
}

# Plot decision probability for 'H0' on y axis and sample size on the x axis
# `summary` is either the data.frame returned by realistic_sim_fixed_size_hill_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_fixed_size_hill_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  pdf("lab-project/figures/report/fixed-decision-prob-delta-compare.pdf", width = 18, height = 6)

  par(mfrow = c(1, 3), mar = c(5, 5, 5, 5))
  mu_vals <- unique(summary_df$mu)
  for (mu in mu_vals) {
    # setup empty plot
    plot(0, 0,
      xlim = c(0, 500), ylim = c(0, 1), type = "n",
      main = bquote(delta == .(mu)),
      ylab = bquote("Decision Probability for " * H[0]),
      xlab = "fixed sample size / stop count n",
      cex.main = 2, cex.lab = 2, cex.axis = 2
    )
    fixed_size <- data.frame(
      trial_count = summary_df$n[summary_df$series == "fixed_size" & summary_df$mu == mu],
      prob = summary_df$prob[summary_df$series == "fixed_size" & summary_df$mu == mu]
    )
    opt_stop <- data.frame(
      stop_count = summary_df$n[summary_df$series == "opt_stop" & summary_df$mu == mu],
      prob = summary_df$prob[summary_df$series == "opt_stop" & summary_df$mu == mu],
      count = summary_df$count[summary_df$series == "opt_stop" & summary_df$mu == mu]
    )
    opt_prob_mean <- sum(opt_stop$prob * opt_stop$count) / 20000

    lines(fixed_size$trial_count, fixed_size$prob, col = my_colors[1], lwd = 5)
    abline(h = opt_prob_mean, col = my_colors[2], lwd = 5, lty = 2)
    # show all lines for mu = 0.02
    if (mu == 0.02) {
      for (mu_all in mu_vals) {
        interp_val <- summary_df$prob[summary_df$series == "stop_avg" & abs(summary_df$mu - mu_all) < 1e-6]
        if (length(interp_val) == 1) {
          abline(h = interp_val, col = my_colors[3 + which(mu_vals == mu_all)], lwd = ifelse(abs(mu_all - 0.02) < 1e-6, 5, 2), lty = c(3 + which(mu_vals == mu_all)))
        }
      }
      legend("bottomright",
        legend = c("fixed n", "optional stopping", bquote("STOP AVG for " * delta == 0.01), bquote("STOP AVG for " * delta == 0.02), bquote("STOP AVG for " * delta == 0.06)),
        col = c(my_colors[1:2], my_colors[4:6]), lwd = 5, lty = c(1, 2, 3 + which(mu_vals == mu_all))
      )
    } else {
      interp_val <- summary_df$prob[summary_df$series == "stop_avg" & abs(summary_df$mu - mu) < 1e-6]
      if (length(interp_val) == 1) {
        abline(h = interp_val, col = my_colors[3 + which(mu_vals == mu)], lwd = 5, lty = c(3 + which(mu_vals == mu)))
      }
    }
  }
  dev.off()
  invisible(summary_df)
}

plot_realistic_sim_fixed_size_hill_summary(realistic_sim_fixed_size_hill_summary(
  fixed_opt_avg_data,
  summary_path = "lab-project/summaries/fixed-decision-prob-delta-compare.txt"
))


###### Prior sampling and weight calculation ######

# Given the Optional Stopping type I and type II error rates,
# we calculate the prior distribution of the effect size.
# The cauchy distribution and the given scale is known.
# This script calculates the weight of the cauchy distribution that is needed to achieve the given error rates.

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

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved prior distribution summary to:", summary_path, "\n")
  }

  summary_df
}

# `summary` is either the data.frame returned by prior_distr_summary()
# or a path to the summary .txt it wrote.
plot_prior_distr_summary <- function(summary) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  bf_crit <- summary_df$bf_crit[1]
  r_val <- summary_df$r_val[1]
  df_type <- summary_df$df_type[1]
  deltas_plot <- summary_df$delta
  probs_plot <- summary_df$prob
  switch(df_type,
    "opt_stop" = {
      col <- my_colors[2]
      title_text <- "Optional Stopping"
    },
    "fixed_max" = {
      col <- my_colors[3]
      title_text <- "MAX"
    },
    "fixed_opt_avg" = {
      col <- my_colors[4]
      title_text <- "STOP AVG"
    },
    "weighted_sum" = {
      col <- my_colors[5]
      title_text <- "SAME DIST"
    },
    stop("Invalid df")
  )

  # 5. Plotting
  pdf(paste0("lab-project/figures/prior-distribution_", df_type, ".pdf"))

  # Set margins to ensure labels aren't cut off
  par(mar = c(5, 5, 4, 2))

  plot(deltas_plot, probs_plot,
    type = "l", col = col, lwd = 3,
    main = bquote("Decision distribution for " * H[0] * " (" * .(title_text) * ")" * ", " * r * "=" * .(round(r_val, 3)) * ", " * BF[crit] * "=" * .(bf_crit)),
    xlab = bquote("Effect size " * delta),
    ylab = "Decision probability",
    cex.main = 0.9
  )

  total_mass <- sum(probs_plot)
  # Dynamic positioning for text so it doesn't overlap the curve
  text(min(deltas_plot) * 0.7, max(probs_plot) * 0.9,
    labels = bquote("Total Mass: " * .(round(total_mass, 3))), cex = 1.2
  )

  dev.off()
  invisible(summary_df)
}

con <- dbConnect(duckdb(), db_file)
plot_prior_distr_summary(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "opt_stop",
  summary_path = "lab-project/summaries/prior-distribution_opt_stop.txt"
))
plot_prior_distr_summary(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "fixed_max",
  summary_path = "lab-project/summaries/prior-distribution_fixed_max.txt"
))
plot_prior_distr_summary(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "fixed_opt_avg",
  summary_path = "lab-project/summaries/prior-distribution_fixed_opt_avg.txt"
))
plot_prior_distr_summary(prior_distr_summary(
  con, bf_crit = 3, r_val = r_vals[2], df_type = "weighted_sum",
  summary_path = "lab-project/summaries/prior-distribution_weighted_sum.txt"
))
dbDisconnect(con)
