# Plots comparing fixed sample size designs with optional stopping (lab project).
# Only reads lab-project/summaries/ (written by lab-project/summarise/fixed-size.R).
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
source("shared/define_colors.R")
source("shared/summary-io.R")
# Initialise params
r_vals <- c(0.5, 1, 2) / sqrt(2)

# Plot decision probability for 'H0' on y axis and sample size on the x axis
# `summary` is either the data.frame returned by realistic_sim_fixed_size_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_fixed_size_summary <- function(summary) {
  summary_df <- read_summary(summary)
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

# `summary` is either the data.frame returned by realistic_sim_fixed_max_summary()
# or a path to the summary .txt it wrote. Draws into the current device.
plot_realistic_sim_fixed_max_summary <- function(summary) {
  summary_df <- read_summary(summary)
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

# `summary` is either the data.frame returned by realistic_sim_combined_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_combined_summary <- function(summary) {
  summary_df <- read_summary(summary)
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

# `summary` is either the data.frame returned by decision_sim_fixed_size_summary()
# or a path to the summary .txt it wrote. Writes one pdf per decision.
plot_decision_sim_fixed_size_summary <- function(summary) {
  summary_df <- read_summary(summary)
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

# Plot decision probability for 'H0' on y axis and sample size on the x axis
# `summary` is either the data.frame returned by realistic_sim_fixed_size_hill_summary()
# or a path to the summary .txt it wrote.
plot_realistic_sim_fixed_size_hill_summary <- function(summary) {
  summary_df <- read_summary(summary)
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

# `summary` is either the data.frame returned by prior_distr_summary()
# or a path to the summary .txt it wrote.
plot_prior_distr_summary <- function(summary) {
  summary_df <- read_summary(summary)
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

plot_realistic_sim_fixed_size_summary("lab-project/summaries/fixed-decision-STOP_AVG.txt")
pdf("lab-project/figures/realistic-sim-fixed-max-bf3.pdf")
plot_realistic_sim_fixed_max_summary("lab-project/summaries/realistic-sim-fixed-max-bf3.txt")
dev.off()
plot_realistic_sim_combined_summary("lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-1.414.txt")
plot_realistic_sim_combined_summary("lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-0.707.txt")
plot_realistic_sim_combined_summary("lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-3-r-0.354.txt")
plot_realistic_sim_combined_summary("lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-6-r-0.707.txt")
plot_realistic_sim_combined_summary("lab-project/summaries/realistic-fixed-sim-overview-for-bfcrit-10-r-0.707.txt")
plot_decision_sim_fixed_size_summary("lab-project/summaries/realistic-sim-fixed-size-all-decisions-bf-crit-3-r-0.707.txt")
plot_realistic_sim_fixed_size_hill_summary("lab-project/summaries/fixed-decision-prob-delta-compare.txt")
plot_prior_distr_summary("lab-project/summaries/prior-distribution_opt_stop.txt")
plot_prior_distr_summary("lab-project/summaries/prior-distribution_fixed_max.txt")
plot_prior_distr_summary("lab-project/summaries/prior-distribution_fixed_opt_avg.txt")
plot_prior_distr_summary("lab-project/summaries/prior-distribution_weighted_sum.txt")
