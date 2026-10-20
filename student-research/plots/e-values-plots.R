# Plots for the e-value simulations of student-research/calculations/e-values-simulation.R.
# Only reads the summaries in student-research/summaries/, so no simulation has to be run.

# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()

alpha <- 0.05

# ==============================================================================
# 3. Plotting Random Walks by Profile
# ==============================================================================
plot_walk_summary <- function(summary_df, eCrit_upper, eCrit_lower,
                              main_title, filename, show_band = TRUE) {
  out_dir <- dirname(filename)
  if (!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  pdf(filename, width = 10, height = 7.5)
  on.exit(dev.off())

  log_upper <- log(eCrit_upper)
  log_lower <- log(eCrit_lower)

  finite_vals <- c(summary_df$log_q25, summary_df$log_q75, log_lower, log_upper)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  if (length(finite_vals) == 0) {
    x_min <- log_lower - 1; x_max <- log_upper + 1
  } else {
    x_min <- min(finite_vals) - 1; x_max <- max(finite_vals) + 1.5
  }

  n_max     <- summary_df$n_max[1]
  trueDelta <- summary_df$trueDelta[1]
  p_opt     <- summary_df$p_optstop_h1[1]
  p_fix     <- summary_df$p_fixed_h1[1]

  # what we condition on: H0 only when the data really come from delta = 0
  cond <- if (isTRUE(all.equal(trueDelta, 0))) "H0" else sprintf("delta=%.2f", trueDelta)

  plot(1, type = "n",
       xlim = c(x_min, x_max),
       ylim = c(0, max(summary_df$n, na.rm = TRUE)),
       xlab = "Log E-Value (0 = Indecisive)",
       ylab = "Number of Observations (n)",
       main = main_title)

  prof_colors <- c(
    "OptStop H1, n_max H1"                 = "forestgreen",
    "OptStop H1, n_max H0"                 = "orange",
    "OptStop H1, n_max indecisive"         = "lightgreen",
    "OptStop H0, n_max H0"                 = "firebrick",
    "OptStop H0, n_max H1"                 = "purple",
    "OptStop H0, n_max indecisive"         = "salmon",
    "OptStop indecisive, n_max H1"         = "steelblue",
    "OptStop indecisive, n_max H0"         = "lightblue",
    "OptStop indecisive, n_max indecisive" = "gray50"
  )

  prof_pct <- tapply(summary_df$profile_pct, summary_df$walk_profile, function(x) x[1])
  prof_pct <- sort(prof_pct, decreasing = TRUE)

  legend_labels <- character(0)
  legend_cols   <- character(0)

  for (prof in names(prof_pct)) {
    df_sub <- summary_df[summary_df$walk_profile == prof, ]
    df_sub <- df_sub[order(df_sub$n), ]
    col_base <- if (prof %in% names(prof_colors)) prof_colors[[prof]] else "black"

    if (show_band) {
      polygon(c(df_sub$log_q25, rev(df_sub$log_q75)),
              c(df_sub$n, rev(df_sub$n)),
              col = adjustcolor(col_base, alpha.f = 0.2), border = NA)
    }
    lines(df_sub$log_median, df_sub$n, col = col_base, lwd = 2)

    # "opt.stop = H1 & fixed-n = H0 : 1.28%"
    parts <- strsplit(sub("^OptStop ", "", prof), ", n_max ")[[1]]
    legend_labels <- c(legend_labels,
                       sprintf("opt.stop = %-10s & fixed-n = %-10s : %5.2f%%",
                               parts[1], parts[2], 100 * prof_pct[[prof]]))
    legend_cols <- c(legend_cols, col_base)
  }

  abline(v = 0, col = "gray50", lty = 3)
  abline(v = log_upper, col = "forestgreen", lty = 2, lwd = 2)
  abline(v = log_lower, col = "firebrick", lty = 2, lwd = 2)

  legend("topleft",
         title = sprintf("Joint decisions: P(opt.stop, fixed-n | %s)", cond),
         legend = c(legend_labels,
                    sprintf("H1 bound for both arms (e = %g)", eCrit_upper),
                    sprintf("H0 bound, opt.stop only (e = %g)", eCrit_lower)),
         col = c(legend_cols, "forestgreen", "firebrick"),
         lty = c(rep(1, length(legend_labels)), 2, 2),
         lwd = c(rep(2, length(legend_labels)), 2, 2),
         bg = "white", cex = 0.7, text.font = 1)

  legend("bottomright",
         legend = c(
           sprintf("P('H1' | %s, opt.stop) = %.4f", cond, p_opt),
           sprintf("P('H1' | %s, fixed-n)  = %.4f", cond, p_fix)
         ),
         bg = "white", cex = 0.8, bty = "o")
  invisible(summary_df)
}
# ==============================================================================
# 5. Plot: sample size (x) vs. decision probability (y)
# ==============================================================================
# `summary` is either the data.frame returned by simulate_fixed_n() or a path
# to the summary .txt it wrote. Draws one line per decision (H1, H0, indecisive).
plot_fixed_n_summary <- function(summary, main_title = NULL, filename = NULL) {
  summary_df <- if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    summary
  }
  summary_df <- summary_df[order(summary_df$n), ]

  if (!is.null(filename)) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    png(filename, width = 900, height = 700, res = 100)
    on.exit(dev.off())
  }

  if (is.null(main_title)) {
    main_title <- sprintf("Decision probabilities by sample size (True Delta = %.2f, %d reps each)",
                          summary_df$trueDelta[1], summary_df$n_reps[1])
  }

  plot(1, type = "n",
       xlim = range(summary_df$n),
       ylim = c(0, 1),
       xlab = "Sample Size (n)",
       ylab = "Decision Probability",
       main = main_title)

  lines(summary_df$n, summary_df$h1_prob, col = "forestgreen", lwd = 2)
  lines(summary_df$n, summary_df$h0_prob, col = "firebrick", lwd = 2)
  lines(summary_df$n, summary_df$ind_prob, col = "gray50", lwd = 2)

  legend("right",
         legend = c(sprintf("H1 (e >= %g)", summary_df$eCrit_upper[1]),
                    sprintf("H0 (e <= %g)", summary_df$eCrit_lower[1]),
                    "Indecisive"),
         col = c("forestgreen", "firebrick", "gray50"),
         lty = 1, lwd = 2,
         bg = "white", cex = 0.8)

  invisible(summary_df)
}

# ==============================================================================
# Simulation 1: Extended Optional Stopping (Variable trueDelta)
# ==============================================================================
eCrit1 <- 1 / alpha
eCrit0 <- alpha
standard_summary <- read.table("student-research/summaries/extended-opt-stop-summary.txt",
                               sep = "\t", header = TRUE, stringsAsFactors = FALSE)
for (trueDelta in unique(standard_summary$trueDelta)) {
  plot_file <- sprintf("student-research/figures/walks-plot-standard-trueDelta-%.2f.pdf", trueDelta)
  plot_walk_summary(
    standard_summary[standard_summary$trueDelta == trueDelta, ],
    eCrit_upper = eCrit1, eCrit_lower = eCrit0,
    main_title = sprintf("Walk Profiles (Alpha bounds | True Delta = %.2f)", trueDelta),
    filename = plot_file
  )
}

# ==============================================================================
# Simulation 2: Asymmetrical Thresholds
# ==============================================================================
asym_upper <- 10        # upper bound, used by BOTH arms
asym_lower <- 0.95    # H0 bound, optional stopping only
asym_summary <- read.table("student-research/summaries/asym-opt-stop-summary.txt",
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
for (n_max in unique(asym_summary$n_max)) {
  sum_df_asym <- asym_summary[asym_summary$n_max == n_max, ]
  plot_file <- sprintf("student-research/figures/walks-plot-asym-n-%.2f.pdf", n_max)
  plot_walk_summary(
    sum_df_asym,
    eCrit_upper = asym_upper, eCrit_lower = asym_lower,
    main_title = sprintf("Walk Profiles (Asymmetric Bounds | True Delta = %.2f)", sum_df_asym$trueDelta[1]),
    filename = plot_file,
    show_band = TRUE
  )
}

# ==============================================================================
# Simulation 3: Fixed sample sizes
# ==============================================================================
for (trueDelta in c("0.00", "0.10", "0.20", "0.30", "0.42", "0.50")) {
  plot_fixed_n_summary(
    sprintf("student-research/summaries/fixed-n-summary-trueDelta-%s-own.txt", trueDelta),
    filename = sprintf("student-research/figures/fixed-n-plot-trueDelta-%s-own.png", trueDelta)
  )
}
