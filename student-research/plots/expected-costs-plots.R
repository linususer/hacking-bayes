# Plot functions for the expected cost analysis. They draw from the summaries
# returned by expected-costs-summaries.R or from the saved summary files.
# Sourced by student-research/plots/expected-costs.R.

# ====================================================================
# Plot function for expected costs (original plot)
# ====================================================================
plot_expected_costs_from_summary <- function(summary, scenario_name = "default") {
    bf_crit <- summary$bf_crit
    r_val <- summary$r_val
    c00 <- summary$c00
    c10 <- summary$c10
    ci0 <- summary$ci0
    c01 <- summary$c01
    c11 <- summary$c11
    ci1 <- summary$ci1
    os_mu <- summary$os_mu
    os_cost <- summary$os_cost
    fm_mu <- summary$fm_mu
    fm_cost <- summary$fm_cost
    fa_delta <- summary$fa_delta
    fa_cost <- summary$fa_cost
    ws_delta <- summary$ws_delta
    ws_cost <- summary$ws_cost

    r_text <- if (abs(r_val - r_vals[1]) < 1e-9) {
        bquote(r == 0.5 / sqrt(2))
    } else if (abs(r_val - r_vals[2]) < 1e-9) {
        bquote(r == 1 / sqrt(2))
    } else {
        bquote(r == 2 / sqrt(2))
    }

    pdf(sprintf(
        "student-research/figures/expected-cost-%s-bf%s-r%.3f.pdf",
        scenario_name, bf_crit, r_val
    ), width = 12)
    par(mar = c(5, 6, 7, 5))

    y_max <- max(c(os_cost, fm_cost, fa_cost, ws_cost), na.rm = TRUE)

    plot(0, 0,
        xlim = c(0, 1), ylim = c(0, 1), type = "n",
        # main = bquote("Expected cost — " * BF[crit] == .(bf_crit) * ",  " * .(r_text)),
        ylab = "Expected cost",
        xlab = bquote("Effect size  " * delta),
        cex.main = 2, cex.lab = 2, cex.axis = 2
    )

    abline(v = 0, col = "grey70", lty = 2, lwd = 1)

    lines(os_mu, os_cost, col = my_colors[2], lwd = 5, lty = 1)
    lines(fm_mu, fm_cost, col = my_colors[3], lwd = 5, lty = 3)
    lines(fa_delta, fa_cost, col = my_colors[4], lwd = 5, lty = 4)
    lines(ws_delta, ws_cost, col = my_colors[5], lwd = 5, lty = 6)

    # Cost-matrix annotation above the plot area
    mtext(
        sprintf(
            "H0 true:  c(H0) = %.2g    c(H1) = %.2g    c(ind) = %.2g",
            c00, c10, ci0
        ),
        side = 3, line = 3.8, cex = 1.1
    )
    mtext(
        sprintf(
            "H1 true:  c(H0) = %.2g    c(H1) = %.2g    c(ind) = %.2g",
            c01, c11, ci1
        ),
        side = 3, line = 2.5, cex = 1.1
    )

    legend("topright",
        legend = c("optional stopping", "MAX", "STOP AVG", "SAME DIST"),
        col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
        lwd = 5, lty = c(1, 3, 4, 6),
        cex = 1.5
    )

    dev.off()
}

# ====================================================================
# Original function that calls the two new ones
# ====================================================================
plot_expected_costs <- function(
    bf_crit,
    r_val,
    cost_matrix = matrix(
        c(
            0, 1, 0,
            1, 0, 0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(
            c("H0_true", "H1_true"),
            c("dec_H0", "dec_H1", "dec_ind")
        )
    ),
    db_file = "data/hacking-bayes.duckdb",
    scenario_name = "default") {
    summary <- summary_expected_costs_original(bf_crit, r_val, cost_matrix, db_file)
    plot_expected_costs_from_summary(summary, scenario_name)
}

plot_cost_by_bf <- function(
    cost_matrix = matrix(
        c(
            0, 1, 0,
            1, 0, 0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(
            c("H0_true", "H1_true"),
            c("dec_H0", "dec_H1", "dec_ind")
        )
    ),
    scenario_name = "default") {
    c00 <- cost_matrix["H0_true", "dec_H0"]
    c10 <- cost_matrix["H0_true", "dec_H1"]
    ci0 <- cost_matrix["H0_true", "dec_ind"]
    c01 <- cost_matrix["H1_true", "dec_H0"]
    c11 <- cost_matrix["H1_true", "dec_H1"]
    ci1 <- cost_matrix["H1_true", "dec_ind"]

    thresholds <- bf_crit_from_cost(cost_matrix)
    bf_crit <- thresholds$bf_crit
    bf_crit2 <- thresholds$bf_crit2

    # BF_01 axis: fine grid on log scale so both thresholds are visible
    bf_seq <- exp(seq(log(0.01), log(100), length.out = 1000))

    C_H0 <- (c00 * bf_seq + c01) / (1 + bf_seq)
    C_H1 <- (c10 * bf_seq + c11) / (1 + bf_seq)
    C_ind <- (ci0 * bf_seq + ci1) / (1 + bf_seq)

    y_max <- max(c(C_H0, C_H1, C_ind), na.rm = TRUE)

    region_col_h1 <- rgb(0.85, 0.90, 1.00, alpha = 0.6) # decide H1  (low BF_01)
    region_col_ind <- rgb(0.93, 0.93, 0.93, alpha = 0.8) # indecisive (between thresholds)
    region_col_h0 <- rgb(1.00, 0.90, 0.85, alpha = 0.6) # decide H0  (high BF_01)
    crit_col <- "black"
    crit_lwd <- 3.5

    x_lo <- 0
    x_hi <- 15 # matches xlim below

    pdf(sprintf("student-research/figures/cost-by-bf-%s.pdf", scenario_name), width = 12)
    par(mar = c(9, 6, 7, 5)) # increase bottom margin from 5 to 9
    plot(0, 0,
        xlim = range(0, 15), ylim = c(0, y_max * 1.1), type = "n",
        # log  = "x",
        main = "",
        ylab = "Expected cost C(d, BF_01)",
        xlab = bquote(BF[`01`]),
        cex.lab = 2, cex.axis = 2
    )

    usr <- par("usr")
    y_bottom <- usr[3]
    y_top <- usr[4]

    # region boundaries clipped to the visible x-range
    bf_crit_clip <- min(max(bf_crit, x_lo), x_hi)
    bf_crit2_clip <- min(max(bf_crit2, x_lo), x_hi)

    rect(x_lo, y_bottom, bf_crit2_clip, y_top, col = region_col_h1, border = NA)
    rect(bf_crit2_clip, y_bottom, bf_crit_clip, y_top, col = region_col_ind, border = NA)
    rect(bf_crit_clip, y_bottom, x_hi, y_top, col = region_col_h0, border = NA)

    # region labels near the top of the plot
    label_y <- y_top - 0.05 * (y_top - y_bottom)
    if (bf_crit2_clip > x_lo + 0.3) {
        text(mean(c(x_lo, bf_crit2_clip)), label_y, "H1",
            cex = 1.4, font = 2, col = "grey20"
        )
    }
    if (bf_crit_clip > bf_crit2_clip + 0.3) {
        text(mean(c(bf_crit2_clip, bf_crit_clip)), label_y, "Ind.",
            cex = 1.4, font = 2, col = "grey20"
        )
    }
    if (x_hi > bf_crit_clip + 0.3) {
        text(mean(c(bf_crit_clip, x_hi)), label_y, "H0",
            cex = 1.4, font = 2, col = "grey20"
        )
    }

    # redraw box + axes on top of the shaded rectangles
    box()

    lines(bf_seq, C_H0, col = my_colors[1], lwd = 5, lty = 1)
    lines(bf_seq, C_H1, col = my_colors[2], lwd = 5, lty = 2)
    lines(bf_seq, C_ind, col = my_colors[3], lwd = 5, lty = 3)

    # ---- Threshold verticals: made more prominent (thicker, solid, dark, labeled) ----
    abline(v = bf_crit, col = crit_col, lty = 1, lwd = crit_lwd)
    abline(v = bf_crit2, col = crit_col, lty = 1, lwd = crit_lwd)

    # inline value labels right next to each threshold line
    # text(bf_crit, y_top - 0.02 * (y_top - y_bottom),
    #   sprintf("BF_crit = %.3f", bf_crit),
    #   pos = 4, cex = 1.2, font = 2, col = crit_col, xpd = NA)
    # text(bf_crit2, y_top - 0.02 * (y_top - y_bottom),
    #   sprintf("BF_crit2 = %.4f", bf_crit2),
    #   pos = 2, cex = 1.2, font = 2, col = crit_col, xpd = NA)

    mtext(
        sprintf(
            "H0 true:  c(H0) = %.2g    c(H1) = %.2g    c(ind) = %.2g",
            c00, c10, ci0
        ),
        side = 3, line = 2.5, cex = 1.1
    )
    mtext(
        sprintf(
            "H1 true:  c(H0) = %.2g    c(H1) = %.2g    c(ind) = %.2g",
            c01, c11, ci1
        ),
        side = 3, line = 1.2, cex = 1.1
    )

    legend(
        x = par("usr")[1], y = par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.18,
        legend = c("C('H_0', BF_01)", "C('H_1', BF_01)", "C('~', BF_01)"),
        col = c(my_colors[1], my_colors[2], my_colors[3]),
        lwd = 5, lty = c(1, 2, 3),
        cex = 1.5,
        horiz = TRUE,
        bty = "n",
        xpd = NA
    )

    dev.off()
}

# `summary` is either the data.frame returned by prior_distr_cost_summary()
# or a path to the summary .txt it wrote.
plot_prior_distr_cost_summary <- function(summary) {
    summary_df <- if (is.character(summary)) {
        read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
        summary
    }
    bf_crit <- summary_df$bf_crit[1]
    r_val <- summary_df$r_val[1]
    df_type <- summary_df$df_type[1]
    scenario_name <- summary_df$scenario_name[1]
    grid_extent_mode <- summary_df$grid_extent_mode[1]
    deltas_plot <- summary_df$delta
    cost_plot <- summary_df$weighted_cost
    h0_point_cost <- summary_df$h0_point_cost[1]
    h1_total_cost <- summary_df$h1_total_cost[1]
    expected_cost <- summary_df$expected_cost[1]

    col <- switch(df_type,
        "opt_stop" = my_colors[2],
        "fixed_max" = my_colors[3],
        "stop_avg" = my_colors[4],
        "same_dist" = my_colors[5],
        stop("Invalid df_type")
    )
    title_text <- switch(df_type,
        "opt_stop" = "Optional Stopping",
        "fixed_max" = "MAX",
        "stop_avg" = "STOP AVG",
        "same_dist" = "SAME DIST",
        stop("invalid df_type")
    )

    # NEW: mode goes into the filename so interval/threshold don't overwrite each other
    pdf(sprintf(
        "student-research/figures/prior-weighted-cost_%s-%s-bf%s-%s.pdf",
        df_type, scenario_name, bf_crit, grid_extent_mode
    ))
    par(mar = c(5, 5, 4, 2))

    y_max <- max(cost_plot, h0_point_cost, na.rm = TRUE)

    plot(deltas_plot, cost_plot,
        type = "l", col = col, lwd = 3,
        ylim = c(0, y_max * 1.15),
        main = bquote("Prior-weighted expected cost (" * .(title_text) * ")" * ", "
            * r * "=" * .(round(r_val, 3)) * ", " * BF[crit] * "=" * .(bf_crit)),
        xlab = bquote("Effect size " * delta),
        ylab = "Prior-weighted expected cost",
        cex.main = 0.9
    )

    # NEW: point-prior cost drawn as a vertical spike at delta = 0
    # (a point mass has no width, so a vertical line is the natural way
    # to represent its "density" alongside the continuous H1 curve)
    segments(
        x0 = 0, y0 = 0, x1 = 0, y1 = h0_point_cost,
        col = "black", lwd = 3, lty = 2
    )
    points(0, h0_point_cost, pch = 19, col = "black")

    legend("topright",
        legend = c("cauchy weighted", "point weighted"),
        col = c(col, "black"), lwd = 3, lty = c(1, 2),
        pch = c(NA, 19), bty = "n", cex = 0.9
    )

    # NEW: three texts as requested
    text(min(deltas_plot) * 0.7, y_max * 0.95,
        labels = bquote("Cost given H0: " * .(round(h0_point_cost, 3))),
        cex = 1.1, adj = 0
    )
    text(min(deltas_plot) * 0.7, y_max * 0.87,
        labels = bquote("Cost given H1: " * .(round(h1_total_cost, 3))),
        cex = 1.1, adj = 0
    )
    text(min(deltas_plot) * 0.7, y_max * 0.79,
        labels = bquote("Expected cost: " * .(round(expected_cost, 3))),
        cex = 1.1, adj = 0
    )

    dev.off()
    invisible(summary_df)
}

# `summary` is the flat table written to expected_costs_bf-<bf_crit>.txt
# (one row per bin), a path to that file, or the list returned by
# summary_expected_costs(). Plots the rows of one df_type and bound.
plot_cost_bound <- function(summary, df_type = NULL, bound = NULL) {
    summary_df <- if (is.character(summary)) {
        read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
        as.data.frame(summary)
    }
    if (is.null(df_type)) df_type <- summary_df$df_type[1]
    if (is.null(bound)) bound <- summary_df$bound[1]
    summary_df <- summary_df[summary_df$df_type == df_type & summary_df$bound == bound, ]

    # Fix color an title text depending on dataframe type
    col <- switch(df_type,
        "opt_stop" = my_colors[2],
        "fixed_max" = my_colors[3],
        "stop_avg" = my_colors[4],
        "same_dist" = my_colors[5],
        stop("Invalid df_type")
    )
    title_text <- switch(df_type,
        "opt_stop" = "Optional Stopping",
        "fixed_max" = "MAX",
        "stop_avg" = "STOP AVG",
        "same_dist" = "SAME DIST",
        stop("Invalid df_type")
    )
    bf_crit <- summary_df$bf_crit[1]
    r_val <- summary_df$r_val[1]
    h0_point_cost <- summary_df$h0_point_cost[1]
    h1_total_cost <- summary_df$h1_total_cost[1]
    expected_cost <- summary_df$expected_cost[1]
    ## Check if really needed
    ## interval_mass <- summary$interval_mass
    ## residual_mass <- summary$residual_mass

    # Get bins, costs and the correct symmetrical order
    # of the bins across all effect sizes
    bins <- data.frame(left = summary_df$bins.left, cost = summary_df$bins.cost)
    step_delta <- bins$left
    step_cost <- bins$cost
    step_delta_full <- c(-rev(step_delta), step_delta)
    step_cost_full <- c(rev(step_cost), step_cost)

    # Plot residuals
    pdf(sprintf(
        "student-research/figures/prior-weighted-cost_%s-%s-bf%s.pdf",
        df_type, bound, bf_crit
    ))
    par(mar = c(5, 5, 4, 2))
    plot(step_delta_full, step_cost_full,
        type = "s", col = col, lwd = 3,
        ylim = c(0, 1),
        main = bquote("Worst-case, Cauchy-weighted cost (" * .(title_text) * ")" * ", "
            * r * "=" * .(round(r_val, 3)) * ", " * BF[crit] * "=" * .(bf_crit)),
        xlab = bquote("Effect size " * delta),
        ylab = "Cost per bin",
        cex.main = 0.9
    )

    segments(
        x0 = 0, y0 = 0, x1 = 0, y1 = h0_point_cost,
        col = "black", lwd = 3, lty = 2
    )
    points(0, h0_point_cost, pch = 19, col = "black")

    legend("topright",
        legend = c(paste(bound, " cauchy costs"), "point costs"),
        col = c(col, "black"), lwd = 3, lty = c(1, 2),
        pch = c(NA, 19), bty = "n", cex = 0.9
    )

    text(min(step_delta_full) * 0.7, 1 * 0.95,
        labels = bquote("Cost given H0: " * .(round(h0_point_cost, 3))), cex = 1.1, adj = 0
    )
    text(min(step_delta_full) * 0.7, 1 * 0.87,
        labels = bquote("Cost given H1: " * .(round(h1_total_cost, 3))), cex = 1.1, adj = 0
    )
    text(min(step_delta_full) * 0.7, 1 * 0.79,
        labels = bquote("Expected cost: " * .(round(expected_cost, 3))), cex = 1.1, adj = 0
    )

    dev.off()
}

plot_overview_expected_costs <- function(summary,
                                          file = "student-research/figures/expected_costs_overview.pdf") {

  if (!is.data.frame(summary)) {
    if (is.list(summary) && !is.null(summary$bf_crit)) {
      summary <- as.data.frame(summary)
    } else if (is.list(summary)) {
      summary <- do.call(rbind, lapply(summary, as.data.frame))
    } else {
      stop("`summary` must be a data.frame or a list of summary_expected_costs() results.")
    }
  }

  key_cols <- c("bf_crit", "r_val", "df_type", "bound",
                "h0_point_cost", "h1_total_cost", "expected_cost")
  missing_cols <- setdiff(key_cols, names(summary))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "summary is missing required column(s): %s\nAvailable columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(names(summary), collapse = ", ")
    ))
  }
  df <- unique(summary[key_cols])

  df_h0 <- unique(df[c("bf_crit", "df_type", "h0_point_cost")])
  names(df_h0)[3] <- "cost"
  df_h0$metric <- "H0"

  df_h1 <- df[c("bf_crit", "df_type", "bound", "h1_total_cost")]
  names(df_h1)[4] <- "cost"
  df_h1$metric <- ifelse(df_h1$bound == "upper", "H1_upper", "H1_lower")
  df_h1$bound <- NULL

  df_exp <- df[c("bf_crit", "df_type", "bound", "expected_cost")]
  names(df_exp)[4] <- "cost"
  df_exp$metric <- ifelse(df_exp$bound == "upper", "Ecost_upper", "Ecost_lower")
  df_exp$bound <- NULL

  plot_df <- rbind(df_h0, df_h1, df_exp)

  df_type_labels <- c(opt_stop  = "Optional Stopping",
                       stop_avg  = "Stop AVG",
                       same_dist = "Same Dist")

  df_types  <- intersect(names(df_type_labels), unique(plot_df$df_type))
  bf_levels <- sort(unique(plot_df$bf_crit))
  metrics   <- c("H0", "H1_upper", "H1_lower", "Ecost_upper", "Ecost_lower")

  pch_map <- c(H0 = 21, H1_upper = 25, H1_lower = 24,
               Ecost_upper = 22, Ecost_lower = 23)   # fillable base-R shapes
  col_map <- c(H0 = "black", H1_upper = my_colors[2], H1_lower = my_colors[2],
               Ecost_upper = my_colors[3], Ecost_lower = my_colors[3])

  plot_df$x <- match(plot_df$df_type, df_types)

  pdf(file, width = 7, height = 5.5)
  on.exit(dev.off(), add = TRUE)

  for (bfc in bf_levels) {
    sub_all <- plot_df[plot_df$bf_crit == bfc, ]

    plot(sub_all$x, sub_all$cost, type = "n", xaxt = "n",
         xlab = "", ylab = "Expected cost",
         main = sprintf("Expected costs overview (BF_crit = %s)", bfc),
         xlim = c(0.5, length(df_types) + 0.5))
    axis(1, at = seq_along(df_types), labels = df_type_labels[df_types])
    grid(nx = NA, ny = NULL, col = "grey85", lty = "dotted")

    for (m in metrics) {
      sub <- sub_all[sub_all$metric == m, ]
      if (nrow(sub) == 0) next
      points(sub$x, sub$cost,
             pch = pch_map[m], bg = col_map[m],
             col = "black", cex = 1.6, lwd = 1.2)
    }

    legend("topright",
           legend = c("H0 cost", "H1 cost upper", "H1 cost lower",
                      "E[cost] upper", "E[cost] lower"),
           pch = pch_map[metrics], pt.bg = col_map[metrics],
           col = "black", bty = "n", cex = 0.8)
  }

  invisible(plot_df)
}
