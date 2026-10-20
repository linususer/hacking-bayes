# Plots for the idealised optional stopping simulation (table bf_decision_threshold).
# Only reads thesis/summaries/ (written by thesis/summarise/optional-stopping.R).
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
source("shared/summary-io.R")

# `summary` is either the data.frame returned by sim_histograms_summary()
# or a path to the summary .txt it wrote. Writes one pdf per mu.
plot_sim_histograms_summary <- function(summary, xlim_max = 60, ylimh0 = 2500, ylimh1 = 200) {
  summary_df <- read_summary(summary)
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

# plot the decision probability for H0 given mu and n_start
# `summary` is either the data.frame returned by decision_prob_curve_summary()
# or a path to the summary .txt it wrote.
# n_start: number of starts for the random walk
plot_decision_prob_curve_summary <- function(summary, n_start) {
  summary_df <- read_summary(summary)
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

plot_sim_histograms_summary("thesis/summaries/asym-both.txt", xlim_max = 60, ylimh0 = 2500, ylimh1 = 2500)
# the committed sym figures use the default H1 axis limit
plot_sim_histograms_summary("thesis/summaries/sym-both.txt", xlim_max = 60, ylimh0 = 2500)
for (case in c("asymmetrical", "symmetrical", "both")) {
  plot_decision_prob_curve_summary(sprintf("thesis/summaries/%s-decision-probability.txt", case), n_start = 1)
}
