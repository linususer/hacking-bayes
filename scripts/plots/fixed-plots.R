# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("scripts/plots/define_colors.R")
# Initialise params
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

db_file <- "data/hacking-bayes.duckdb"

####################
##### QUERIES ######
####################
get_opt_stop <- function(con, bf_crit, r_val, decision = 0) {
  dbExecute(con, "SET max_expression_depth TO 10000")
  opt_stop <- dbGetQuery(
    con, "
    SELECT mu, COUNT (CASE WHEN decision = ? THEN 1.0 END) * 1.0 / COUNT(*) AS prob,
           AVG(stop_count) AS mean_count
    FROM cauchy_sym
    WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
      AND trial_start = 2 AND trial_end = 100000
    GROUP BY mu
    ORDER BY mu",
    list(decision, r_val, bf_crit)
  )
  opt_stop
}

get_fixed_max <- function(con, bf_crit, r_val, decision = 0) {
  dbExecute(con, "SET max_expression_depth TO 10000")
  fixed_max <- dbGetQuery(
    con, "
    WITH prob_table AS (
      SELECT mu, trial_count,
             COUNT(CASE WHEN decision = ? THEN 1 END) * 1.0 / COUNT(*) AS prob
      FROM cauchy_sym_fixed_size
      WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
      GROUP BY mu, trial_count
    ),
    ranked AS (
      SELECT mu, trial_count, prob,
             ROW_NUMBER() OVER (
               PARTITION BY mu
               ORDER BY prob DESC, trial_count ASC
             ) AS rn
      FROM prob_table
    )
    SELECT mu, trial_count, prob
    FROM ranked
    WHERE rn = 1
    ORDER BY mu",
    list(decision, r_val, bf_crit)
  )
  fixed_max
}

get_prob_by_trial <- function(con, bf_crit, r_val, decision = 1, fixed_max_df) {
  dbExecute(con, "SET max_expression_depth TO 10000")
  values_list <- paste(
    apply(fixed_max_df, 1, function(row) {
      sprintf("SELECT %s AS mu, %s AS trial_count", row["mu"], row["trial_count"])
    }),
    collapse = "\nUNION ALL\n"
  )

  sql <- sprintf("
    WITH selected_combinations AS (
      %s
    )
    SELECT s.mu, s.trial_count,
           COUNT(CASE WHEN c.decision = ? THEN 1 END) * 1.0 / COUNT(*) AS prob
    FROM cauchy_sym_fixed_size c
    JOIN selected_combinations s
      ON ABS(c.mu - s.mu) < 1e-6 AND c.trial_count = s.trial_count
    WHERE ABS(? - c.r) < 1e-6 AND c.bf_crit = ?
    GROUP BY s.mu, s.trial_count
    ORDER BY s.mu, s.trial_count
  ", values_list)

  params <- list(decision, r_val, bf_crit)

  prob_by_trial <- dbGetQuery(con, sql, params)
  prob_by_trial
}



get_fixed_opt_avg <- function(con, bf_crit, r_val, decision = 0) {
  dbExecute(con, "SET max_expression_depth TO 10000")
  fixed_opt_avg <- dbGetQuery(con, "
    WITH opt AS (
      SELECT mu, AVG(stop_count) AS avg_stop,
             SUM(CASE WHEN decision = ? THEN 1 END)*1.0/COUNT(*) AS prob_opt
      FROM cauchy_sym
      WHERE ABS(r - ?) < 1e-6 AND bf_crit = ? AND trial_start = 2 AND trial_end = 100000
      GROUP BY mu
    ),
    fixed AS (
      SELECT mu, trial_count,
             SUM(CASE WHEN decision = ? THEN 1 ELSE 0 END)*1.0/COUNT(*) AS prob_fixed
      FROM cauchy_sym_fixed_size
      WHERE ABS(r - ?) < 1e-6 AND bf_crit = ?
      GROUP BY mu, trial_count
    ),
    interp AS (
      SELECT
        o.mu AS delta,
        o.avg_stop,
        o.prob_opt,
        f1.trial_count AS t1,
        f1.prob_fixed AS p1,
        f2.trial_count AS t2,
        f2.prob_fixed AS p2,
        CASE
          WHEN f1.trial_count = f2.trial_count THEN f1.prob_fixed
          ELSE f1.prob_fixed + (o.avg_stop - f1.trial_count) * (f2.prob_fixed - f1.prob_fixed) / (f2.trial_count - f1.trial_count)
        END AS prob
      FROM opt o
      JOIN fixed f1 ON ABS(f1.mu - o.mu) < 0.001 AND f1.trial_count <= o.avg_stop
      JOIN fixed f2 ON ABS(f2.mu - o.mu) < 0.001 AND f2.trial_count >= o.avg_stop
      WHERE f1.trial_count = (
        SELECT MAX(trial_count)
        FROM fixed
        WHERE ABS(mu - o.mu) < 0.001 AND trial_count <= o.avg_stop
      )
      AND f2.trial_count = (
        SELECT MIN(trial_count)
        FROM fixed
        WHERE ABS(mu - o.mu) < 0.001 AND trial_count >= o.avg_stop
      )
    )
    SELECT delta, avg_stop, prob_opt, prob
    FROM interp
    ORDER BY delta",
    params = list(decision, r_val, bf_crit, decision, r_val, bf_crit)
  )
  fixed_opt_avg
}

get_fixed_weighted_sum <- function(con, bf_crit, r_val, decision = 0) {
  dbExecute(con, "SET max_expression_depth TO 10000")
  weighted_sum <- dbGetQuery(con, "
WITH stop_dist AS (
  SELECT
    ROUND(mu, 6) AS mu,
    stop_count,
    COUNT(*)::DOUBLE PRECISION AS total
  FROM cauchy_sym
  WHERE ABS(r - ?) < 1e-6
    AND bf_crit = ?
    AND trial_start = 2
    AND trial_end = 100000
  GROUP BY ROUND(mu, 6), stop_count
),
sum_totals AS (
  SELECT
    mu,
    SUM(total) AS total_sum
  FROM stop_dist
  GROUP BY mu
),
fixed_probs AS (
  SELECT
    ROUND(mu, 6) AS mu,
    trial_count,
    AVG(CASE WHEN decision = ? THEN 1.0 ELSE 0.0 END) AS p_h0_fixed
  FROM cauchy_sym_fixed_size
  WHERE ABS(r - ?) < 1e-6
    AND bf_crit = ?
  GROUP BY ROUND(mu, 6), trial_count
),
interpolated AS (
  SELECT
    s.mu,
    s.stop_count,
    s.total,
    st.total_sum,
    -- Find bounding trial_counts
    f1.trial_count AS t1,
    f1.p_h0_fixed AS p1,
    f2.trial_count AS t2,
    f2.p_h0_fixed AS p2,
    CASE
      WHEN f1.trial_count = f2.trial_count THEN f1.p_h0_fixed
      ELSE f1.p_h0_fixed + (s.stop_count - f1.trial_count) * (f2.p_h0_fixed - f1.p_h0_fixed) / (f2.trial_count - f1.trial_count)
    END AS p_interp
  FROM stop_dist s
  JOIN sum_totals st ON s.mu = st.mu
  LEFT JOIN fixed_probs f1 ON f1.mu = s.mu AND f1.trial_count = (
    SELECT MAX(trial_count) FROM fixed_probs
    WHERE mu = s.mu AND trial_count <= s.stop_count
  )
  LEFT JOIN fixed_probs f2 ON f2.mu = s.mu AND f2.trial_count = (
    SELECT MIN(trial_count) FROM fixed_probs
    WHERE mu = s.mu AND trial_count >= s.stop_count
  )
),
weighted AS (
  SELECT
    mu,
    (total / total_sum) * p_interp AS weight_component
  FROM interpolated
  WHERE p_interp IS NOT NULL
)
SELECT
  mu AS delta,
  SUM(weight_component) AS prob
FROM weighted
GROUP BY mu
ORDER BY delta;
", list(r_val, bf_crit, decision, r_val, bf_crit))
  weighted_sum
}

###############
#### PLOTS ####
################

fixed_big_mus <- c(seq(0, 0.09, 0.01), seq(0.1, 1, 0.05))
BF_crits <- c(3)

con <- dbConnect(duckdb(), db_file)
fixed_opt_avg_data <- get_fixed_opt_avg(con, BF_crits[1], r_vals[2])
fixed_weighted_sum_data <- get_fixed_weighted_sum(con, BF_crits[1], r_vals[2])
fixed_max_data <- get_fixed_max(con, BF_crits[1], r_vals[2])
dbDisconnect(con)

# Plot decision probability for 'H0' on y axis and sample size on the x axis
realistic_sim_fixed_size_plot <- function(fixed_opt_avg_data, fixed_weighted_sum_data, fixed_max_data) {
  # increase the size, including font size
  con <- dbConnect(duckdb(), db_file)
  for (mu in c(0.5)) { # c(0.01,0.02,0.06,0.05, 0.1,0.2,0.5)
    pdf(paste("figures/report/fixed-decision-STOP_AVG.pdf", sep = ""))
    par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
    plot(0, 0,
      xlim = c(0, 100), ylim = c(0, 1), type = "n",
      main = bquote("STOP AVG, " * delta == .(mu)),
      ylab = bquote("Decision Probability for " * H[0]), xlab = "fixed sample size / stop count n",
      cex.main = 2, cex.lab = 2, cex.axis = 2
    )
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
      interp_val <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu) < 1e-6, "prob"]
      trial_at_interp_val <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu) < 1e-6, "avg_stop"]
      weighted_val <- fixed_weighted_sum_data[abs(fixed_weighted_sum_data$delta - mu) < 1e-6, "prob"]
      max_val <- fixed_max_data[abs(fixed_max_data$mu - mu) < 1e-6, "prob"]
      trial_at_max_val <- fixed_max_data[abs(fixed_max_data$mu - mu) < 1e-6, "trial_count"]

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
  dbDisconnect(con)
}
realistic_sim_fixed_size_plot(fixed_opt_avg_data, fixed_weighted_sum_data, fixed_max_data)

# Plot maximum probability of 'H0' for each mu

realistic_sim_fixed_max <- function(bf_crit, r) {
  con <- dbConnect(duckdb(), db_file)
  # Decision probability for H0 given mu
  # pdf(paste("figures/realistic-sim-fixed-max-bf3.pdf", sep = ""))
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability for 'H0'"),
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote("Effect size " * delta)
  )
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
  # get the first index where the probability is smaller than 0.5
  lines(fixed_max[["mu"]], fixed_max[["prob"]], col = "black", lwd = 2)
  lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd = 2)
  legend("topright", legend = c("fixed probability maximum", "optional stopping probability"), fill = c("black", my_colors[2]))
  # dev.off()
  dbDisconnect(con)
  max_arg <- max_arg <- merge(
    fixed_max[, c("mu", "trial_count")],
    opt_stop[, c("mu", "mean_count")],
    by = "mu",
    suffixes = c("_fixed", "_opt_stop")
  )
  print(max_arg)
}


# P('Diff') = P(H0 | n_opt_stop_avg) and P(H0 | n_fixed) for n_fixed ~= n_opt_stop_avg
# get stop_counts from optional stopping results
#  P('Diff')


realistic_sim_combined_plot <- function(bf_crit, r_val, db_file = "data/hacking-bayes.duckdb", decision = 0) {
  con <- dbConnect(duckdb(), db_file)

  opt_stop <- get_opt_stop(con, bf_crit, r_val, decision)
  fixed_max <- get_fixed_max(con, bf_crit, r_val, decision)
  fixed_opt_avg <- get_fixed_opt_avg(con, bf_crit, r_val, decision)
  weighted_results <- get_fixed_weighted_sum(con, bf_crit, r_val, decision)

  dbDisconnect(con)

  # Plot
  par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
  pdf(paste("figures/report/realistic-fixed-sim-overview-for-bfcrit-", bf_crit, "-r-", round(r_val, 3), ".pdf", sep = ""))
  par(mar = c(5, 6, 5, 5))
  # r_val text
  if (r_val == r_vals[1]) {
    r_text <- bquote(r == 0.5 / sqrt(2))
  } else if (r_val == r_vals[2]) {
    r_text <- bquote(r == 1 / sqrt(2))
  } else if (r_val == r_vals[3]) {
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
}

realistic_sim_combined_plot(3, r_vals[3])
realistic_sim_combined_plot(3, r_vals[2])
realistic_sim_combined_plot(3, r_vals[1])
realistic_sim_combined_plot(6, r_vals[2])
realistic_sim_combined_plot(10, r_vals[2])

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
#   pdf(paste("figures/stop-count-to-effect-size-r-", round(r_val, 3), ".pdf", sep = ""), width = 12, height = 8)
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
#   pdf(paste("figures/stop-count-to-effect-size-bf-", bf_crit, ".pdf", sep = ""), width = 12, height = 8)
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
#   pdf(paste("figures/fixed-vs-optional-stopping-bf-", bf_crit, "-r-", round(r_val, 3), ".pdf", sep = ""), width = 12, height = 8)
#   plot(0, 0,
#     xlim = c(0, 1), ylim = c(0, 1), type = "n",
#     main = bquote("Fixed Size vs Optional Stopping for " * BF[crit] * " = " * .(bf_crit) * ", r = " * .(r_val)),
#     ylab = bquote("Decision Probability for " * H[0]),
#     xlab = bquote("Effect Size " * delta)
#   )

#   fixed_max <- get_fixed_max(con, bf_crit, 2 * r_val)
#   fixed_opt_avg <- get_fixed_opt_avg(con, bf_crit, 2 * r_val)
#   weighted_sum <- get_fixed_weighted_sum(con, bf_crit, 2 * r_val)
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

plot_decision_sim_fixed_size <- function(bf_crit, r_val) {
  con <- dbConnect(duckdb(), db_file)
  fixed_max_d0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)
  for (decision in 0:2) {
    pdf_filename <- sprintf("figures/report/realistic-sim-fixed-size-all-decisions-bf-crit-%s-r-%.3f-%s-wo-legend.pdf", bf_crit, r_val, decision)
    pdf(pdf_filename)

    decision_text <- switch(as.character(decision),
      "0" = bquote(H[0]),
      "1" = bquote(H[1]),
      "2" = bquote("Indecisive")
    )

    if (decision == 0) {
      fixed_max <- fixed_max_d0
    } else {
      fixed_max <- get_prob_by_trial(con, bf_crit, r_val, decision, fixed_max_df = fixed_max_d0)
    }
    print(fixed_max)
    opt_stop <- get_opt_stop(con, bf_crit, r_val, decision)
    fixed_opt_avg <- get_fixed_opt_avg(con, bf_crit, r_val, decision)
    weighted_results <- get_fixed_weighted_sum(con, bf_crit, r_val, decision)
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
  dbDisconnect(con)
}

plot_decision_sim_fixed_size(3, r_vals[2])

BF_crits <- c(10)
# Plot decision probability for 'H0' on y axis and sample size on the x axis
realistic_sim_fixed_size_hill <- function(fixed_opt_avg_data) {
  con <- dbConnect(duckdb(), db_file)
  pdf("figures/report/fixed-decision-prob-delta-compare.pdf", width = 18, height = 6)

  par(mfrow = c(1, 3), mar = c(5, 5, 5, 5))
  mu_vals <- c(0.01, 0.02, 0.06)
  BF_crit <- 10
  for (mu in mu_vals) {
    # setup empty plot
    plot(0, 0,
      xlim = c(0, 500), ylim = c(0, 1), type = "n",
      main = bquote(delta == .(mu)),
      ylab = bquote("Decision Probability for " * H[0]),
      xlab = "fixed sample size / stop count n",
      cex.main = 2, cex.lab = 2, cex.axis = 2
    )
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
    opt_prob_mean <- sum(opt_stop$prob * opt_stop$count) / 20000

    lines(fixed_size$trial_count, fixed_size$prob, col = my_colors[1], lwd = 5)
    abline(h = opt_prob_mean, col = my_colors[2], lwd = 5, lty = 2)
    # show all lines for mu = 0.02
    if (mu == 0.02) {
      for (mu_all in mu_vals) {
        interp_val <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu_all) < 1e-6, "prob"]
        if (length(interp_val) == 1) {
          abline(h = interp_val, col = my_colors[3 + which(mu_vals == mu_all)], lwd = ifelse(abs(mu_all - 0.02) < 1e-6, 5, 2), lty = c(3 + which(mu_vals == mu_all)))
        }
      }
      legend("bottomright",
        legend = c("fixed n", "optional stopping", bquote("STOP AVG for " * delta == 0.01), bquote("STOP AVG for " * delta == 0.02), bquote("STOP AVG for " * delta == 0.06)),
        col = c(my_colors[1:2], my_colors[4:6]), lwd = 5, lty = c(1, 2, 3 + which(mu_vals == mu_all))
      )
    } else {
      interp_val <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu) < 1e-6, "prob"]
      if (length(interp_val) == 1) {
        abline(h = interp_val, col = my_colors[3 + which(mu_vals == mu)], lwd = 5, lty = c(3 + which(mu_vals == mu)))
      }
    }
  }
  dev.off()
  dbDisconnect(con)
}

realistic_sim_fixed_size_hill(fixed_opt_avg_data)


###### Prior sampling and weight calculation ######

# Given the Optional Stopping type I and type II error rates,
# we calculate the prior distribution of the effect size.
# The cauchy distribution and the given scale is known.
# This script calculates the weight of the cauchy distribution that is needed to achieve the given error rates.

plot_prior_distr <- function(con, bf_crit, r_val, df_type = "opt_stop") {
  switch(df_type,
    "opt_stop" = {
      df <- get_opt_stop(con, bf_crit, r_val, decision = 0)
      col <- my_colors[2]
      target_delta <- 0.01
      delta_or_mu <- "mu"
      title_text <- "Optional Stopping"
    },
    "fixed_max" = {
      df <- get_fixed_max(con, bf_crit, r_val, decision = 0)
      col <- my_colors[3]
      target_delta <- 0.05
      delta_or_mu <- "mu"
      title_text <- "MAX"
    },
    "fixed_opt_avg" = {
      df <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 0)
      col <- my_colors[4]
      target_delta <- 0.05
      delta_or_mu <- "delta"
      title_text <- "STOP AVG"
    },
    "weighted_sum" = {
      df <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 0)
      col <- my_colors[5]
      target_delta <- 0.05
      delta_or_mu <- "delta"
      title_text <- "SAME DIST"
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

  # 5. Plotting
  pdf(paste0("figures/prior-distribution_", df_type, ".pdf"))

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
}

con <- dbConnect(duckdb(), db_file)
plot_prior_distr(con, bf_crit = 3, r_val = r_vals[2], df_type = "opt_stop")
plot_prior_distr(con, bf_crit = 3, r_val = r_vals[2], df_type = "fixed_max")
plot_prior_distr(con, bf_crit = 3, r_val = r_vals[2], df_type = "fixed_opt_avg")
plot_prior_distr(con, bf_crit = 3, r_val = r_vals[2], df_type = "weighted_sum")
dbDisconnect(con)