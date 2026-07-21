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
      AND trial_start = 20 AND trial_end = 200
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
      WHERE ABS(r - ?) < 1e-6 AND bf_crit = ? AND trial_start = 20 AND trial_end = 200
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



#####################################
####### Expected cost plots #########
#####################################

# cost_matrix is a 2x3 list
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
  c00 <- cost_matrix["H0_true", "dec_H0"] # correct H_0
  c10 <- cost_matrix["H0_true", "dec_H1"] # Type I
  ci0 <- cost_matrix["H0_true", "dec_ind"] # indecisive | H_0 true
  c01 <- cost_matrix["H1_true", "dec_H0"] # Type II
  c11 <- cost_matrix["H1_true", "dec_H1"] # correct H_1
  ci1 <- cost_matrix["H1_true", "dec_ind"] # indecisive | H_1 true

  # expected cost function calculation
  # P('H_0') = p0
  # P('H_1') = p1
  # P(indecisive') = p2
  ec <- function(p0, p1, p2, delta_vals) {
    is_h0 <- abs(delta_vals) < 1e-9
    ifelse(is_h0,
      c00 * p0 + c10 * p1 + ci0 * p2, # H_0 true
      c01 * p0 + c11 * p1 + ci1 * p2 # H_1 true
    )
  }
  reagrid <- function(src_x, src_y, target_x) {
    approx(src_x, src_y, xout = target_x, rule = 2)$y
  }


  con <- dbConnect(duckdb(), db_file)

  # Optional stopping with mu / delta as the x-axis
  os0 <- get_opt_stop(con, bf_crit, r_val, decision = 0)
  os1 <- get_opt_stop(con, bf_crit, r_val, decision = 1)
  os2 <- get_opt_stop(con, bf_crit, r_val, decision = 2)
  os_mu <- os0$mu
  os_cost <- ec(
    os0$prob,
    reagrid(os1$mu, os1$prob, os_mu),
    reagrid(os2$mu, os2$prob, os_mu),
    os_mu
  )
  

  # Fixed MAX — decision 0 sets the (mu, trial_count) pairs; 1 & 2 follow
  fm0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)
  fm1 <- get_prob_by_trial(con, bf_crit, r_val, decision = 1, fixed_max_df = fm0)
  fm2 <- get_prob_by_trial(con, bf_crit, r_val, decision = 2, fixed_max_df = fm0)
  fm_mu <- fm0$mu
  fm_cost <- ec(
    fm0$prob,
    reagrid(fm1$mu, fm1$prob, fm_mu),
    reagrid(fm2$mu, fm2$prob, fm_mu),
    fm_mu
  )

  # Fixed STOP AVG — uses delta as the x-axis
  fa0 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 0)
  fa1 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 1)
  fa2 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 2)
  fa_delta <- fa0$delta
  fa_cost <- ec(
    fa0$prob,
    reagrid(fa1$delta, fa1$prob, fa_delta),
    reagrid(fa2$delta, fa2$prob, fa_delta),
    fa_delta
  )

  # Fixed SAME DIST — uses delta as the x-axis
  ws0 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 0)
  ws1 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 1)
  ws2 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 2)
  ws_delta <- ws0$delta
  ws_cost <- ec(
    ws0$prob,
    reagrid(ws1$delta, ws1$prob, ws_delta),
    reagrid(ws2$delta, ws2$prob, ws_delta),
    ws_delta
  )

  dbDisconnect(con)

  # Plot
  r_text <- if (abs(r_val - r_vals[1]) < 1e-9) {
    bquote(r == 0.5 / sqrt(2))
  } else if (abs(r_val - r_vals[2]) < 1e-9) {
    bquote(r == 1 / sqrt(2))
  } else {
    bquote(r == 2 / sqrt(2))
  }

  pdf(sprintf(
    "figures/expected-cost-%s-bf%s-r%.3f.pdf",
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
# Calculate BF_crit from costs
bf_crit_from_cost <- function(cost_matrix) {
  c00 <- cost_matrix["H0_true", "dec_H0"]
  c10 <- cost_matrix["H0_true", "dec_H1"]
  ci0 <- cost_matrix["H0_true", "dec_ind"]
  c01 <- cost_matrix["H1_true", "dec_H0"]
  c11 <- cost_matrix["H1_true", "dec_H1"]
  ci1 <- cost_matrix["H1_true", "dec_ind"]

  # C('H_0', BF01) = C('~', BF01)  ->  upper threshold, decide H_0
  bf_crit <- (ci1 - c01) / (c00 - ci0)

  # C('H_1', BF01) = C('~', BF01)  ->  lower threshold, decide H_1
  bf_crit2 <- (ci1 - c11) / (c10 - ci0)

  if (bf_crit <= 1) {
    warning("BF_crit <= 1: indecision region on H_0 side vanishes")
  }
  if (bf_crit2 >= 1) {
    warning("BF_crit2 >= 1: indecision region on H_1 side vanishes")
  }
  if (abs(bf_crit2 - 1 / bf_crit) > 1e-9) {
    message(sprintf(
      "Asymmetric matrix: BF_crit = %.4f, 1/BF_crit2 = %.4f",
      bf_crit, 1 / bf_crit2
    ))
  }

  list(bf_crit = bf_crit, bf_crit2 = bf_crit2)
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

  # ---- NEW: colors for the three decision-region backgrounds ----
  region_col_h1  <- rgb(0.85, 0.90, 1.00, alpha = 0.6)  # decide H1  (low BF_01)
  region_col_ind <- rgb(0.93, 0.93, 0.93, alpha = 0.8)  # indecisive (between thresholds)
  region_col_h0  <- rgb(1.00, 0.90, 0.85, alpha = 0.6)  # decide H0  (high BF_01)
  # more prominent threshold-line styling
  crit_col <- "black"
  crit_lwd <- 3.5

  x_lo <- 0
  x_hi <- 15 # matches xlim below

  pdf(sprintf("figures/cost-by-bf-%s.pdf", scenario_name), width = 12)
  par(mar = c(9, 6, 7, 5)) # increase bottom margin from 5 to 9
  plot(0, 0,
    xlim = range(0, 15), ylim = c(0, y_max * 1.1), type = "n",
    # log  = "x",
    main = "",
    ylab = "Expected cost C(d, BF_01)",
    xlab = bquote(BF[`01`]),
    cex.lab = 2, cex.axis = 2
  )

  # ---- NEW: shade the three decision regions first, behind everything ----
  usr <- par("usr")
  y_bottom <- usr[3]
  y_top <- usr[4]

  # region boundaries clipped to the visible x-range
  bf_crit_clip  <- min(max(bf_crit, x_lo), x_hi)
  bf_crit2_clip <- min(max(bf_crit2, x_lo), x_hi)

  rect(x_lo, y_bottom, bf_crit2_clip, y_top, col = region_col_h1, border = NA)
  rect(bf_crit2_clip, y_bottom, bf_crit_clip, y_top, col = region_col_ind, border = NA)
  rect(bf_crit_clip, y_bottom, x_hi, y_top, col = region_col_h0, border = NA)

  # region labels near the top of the plot
  label_y <- y_top - 0.05 * (y_top - y_bottom)
  if (bf_crit2_clip > x_lo + 0.3) {
    text(mean(c(x_lo, bf_crit2_clip)), label_y, "H1",
      cex = 1.4, font = 2, col = "grey20")
  }
  if (bf_crit_clip > bf_crit2_clip + 0.3) {
    text(mean(c(bf_crit2_clip, bf_crit_clip)), label_y, "Ind.",
      cex = 1.4, font = 2, col = "grey20")
  }
  if (x_hi > bf_crit_clip + 0.3) {
    text(mean(c(bf_crit_clip, x_hi)), label_y, "H0",
      cex = 1.4, font = 2, col = "grey20")
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
#########################
######## CONFIG #########
#########################

# Cost matrix 1: symmetrical costs and indecisive outcomes cost > Type I & Type II error cost
cost_matrix1 <- matrix(
  c(
    0, 0.5, 1,
    0.5, 0, 1
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("H0_true", "H1_true"),
    c("dec_H0", "dec_H1", "dec_ind")
  )
)
# Cost matrix 2: symmetrical costs and indecisive outcomes = Type 1 & Type II error cost
cost_matrix2 <- matrix(
  c(
    0, 1, 1,
    1, 0, 1
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("H0_true", "H1_true"),
    c("dec_H0", "dec_H1", "dec_ind")
  )
)
# Cost matrix 3: symmetrical costs and indecisive outcomes cost < Type I & Type II error cost
cost_matrix3 <- matrix(
  c(
    0, 1, 0.5,
    1, 0, 0.5
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("H0_true", "H1_true"),
    c("dec_H0", "dec_H1", "dec_ind")
  )
)

# Plots for expected costs in general
plot_expected_costs(3, r_vals[2], scenario_name = "default")
plot_expected_costs(3, r_vals[2], cost_matrix1, scenario_name = "sym-ind-ge-error")
plot_expected_costs(3, r_vals[2], cost_matrix2, scenario_name = "sym-ind-e-error")
plot_expected_costs(3, r_vals[2], cost_matrix3, scenario_name = "sym-ind-le-error")
# Plots for BF_crit
plot_cost_by_bf(scenario_name = "default")
plot_cost_by_bf(cost_matrix1, scenario_name = "sym-ind-ge-error")
plot_cost_by_bf(cost_matrix2, scenario_name = "sym-ind-e-error")
plot_cost_by_bf(cost_matrix3, scenario_name = "sym-ind-le-error")


# Cost matrix 4:
cost_matrix4 <- matrix(
  c(
    0, 1, 0.25,
    1, 0, 0.25
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("H0_true", "H1_true"),
    c("dec_H0", "dec_H1", "dec_ind")
  )
)

plot_expected_costs(3, r_vals[2], cost_matrix4, scenario_name = "sym-ind-bfcrit3")
plot_cost_by_bf(cost_matrix4, scenario_name = "sym-ind-bfcrit3")

cost_matrix5 <- matrix(
  c(
    0,1,1/7,
    1,0,1/7
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("H0_true", "H1_true"),
    c("dec_H0", "dec_H1", "dec_ind")
  )
)

cost_matrix6 <- matrix(
  c(
    0,1,1/11,
    1,0,1/11
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    c("H0_true", "H1_true"),
    c("dec_H0", "dec_H1", "dec_ind")
  )
)

plot_expected_costs(6, r_vals[2], cost_matrix5, scenario_name = "sym-ind-bfcrit6")
plot_cost_by_bf(cost_matrix5, scenario_name = "sym-ind-bfcrit6")
plot_expected_costs(10, r_vals[2], cost_matrix6, scenario_name = "sym-ind-bfcrit10")
plot_cost_by_bf(cost_matrix6, scenario_name = "sym-ind-bfcrit10")

###########################################
############ WEIGHT COSTS #################
###########################################
get_cost_curve <- function(con, bf_crit, r_val, cost_matrix, df_type = "opt_stop") {
  c00 <- cost_matrix["H0_true", "dec_H0"]
  c10 <- cost_matrix["H0_true", "dec_H1"]
  ci0 <- cost_matrix["H0_true", "dec_ind"]
  c01 <- cost_matrix["H1_true", "dec_H0"]
  c11 <- cost_matrix["H1_true", "dec_H1"]
  ci1 <- cost_matrix["H1_true", "dec_ind"]

  ec <- function(p0, p1, p2, delta_vals) {
    is_zero <- abs(delta_vals) < 1e-9   # FIXED: was `is_h0`, undefined at point of use
    ifelse(is_zero,
      c01 * p0 + ci1 * p2,             # delta = 0: H1 ground truth; only dec-H0 and indecisive costs
      c01 * p0 + c11 * p1 + ci1 * p2   # delta != 0: unchanged, full H1-true blend
    )
  }
  reagrid <- function(src_x, src_y, target_x) {
    approx(src_x, src_y, xout = target_x, rule = 2)$y
  }

  switch(df_type,
    "opt_stop" = {
      d0 <- get_opt_stop(con, bf_crit, r_val, decision = 0)
      d1 <- get_opt_stop(con, bf_crit, r_val, decision = 1)
      d2 <- get_opt_stop(con, bf_crit, r_val, decision = 2)
      x <- d0$mu
      data.frame(mu = x, cost = ec(d0$prob, reagrid(d1$mu, d1$prob, x), reagrid(d2$mu, d2$prob, x), x))
    },
    "fixed_max" = {
      d0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)
      d1 <- get_prob_by_trial(con, bf_crit, r_val, decision = 1, fixed_max_df = d0)
      d2 <- get_prob_by_trial(con, bf_crit, r_val, decision = 2, fixed_max_df = d0)
      x <- d0$mu
      data.frame(mu = x, cost = ec(d0$prob, reagrid(d1$mu, d1$prob, x), reagrid(d2$mu, d2$prob, x), x))
    },
    "fixed_opt_avg" = {
      d0 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 0)
      d1 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 1)
      d2 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 2)
      x <- d0$delta
      data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$delta, d1$prob, x), reagrid(d2$delta, d2$prob, x), x))
    },
    "weighted_sum" = {
      d0 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 0)
      d1 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 1)
      d2 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 2)
      x <- d0$delta
      data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$delta, d1$prob, x), reagrid(d2$delta, d2$prob, x), x))
    },
    stop("Invalid df_type")
  )
}

# NEW: cost under the point-null prior (H0 exactly true, weight = 1),
# using the H0-true row of the cost matrix, evaluated at delta = 0.
get_h0_point_cost <- function(con, bf_crit, r_val, cost_matrix, df_type = "opt_stop") {
  c00 <- cost_matrix["H0_true", "dec_H0"]
  c10 <- cost_matrix["H0_true", "dec_H1"]
  ci0 <- cost_matrix["H0_true", "dec_ind"]

  p_at_zero <- function(df, x_col) {
    idx <- which(abs(df[[x_col]]) < 1e-9)
    if (length(idx) == 0) return(0)
    df$prob[idx[1]]
  }

  switch(df_type,
    "opt_stop" = {
      d0 <- get_opt_stop(con, bf_crit, r_val, decision = 0)
      d1 <- get_opt_stop(con, bf_crit, r_val, decision = 1)
      d2 <- get_opt_stop(con, bf_crit, r_val, decision = 2)
      p0 <- p_at_zero(d0, "mu"); p1 <- p_at_zero(d1, "mu"); p2 <- p_at_zero(d2, "mu")
    },
    "fixed_max" = {
      d0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)
      d1 <- get_prob_by_trial(con, bf_crit, r_val, decision = 1, fixed_max_df = d0)
      d2 <- get_prob_by_trial(con, bf_crit, r_val, decision = 2, fixed_max_df = d0)
      p0 <- p_at_zero(d0, "mu"); p1 <- p_at_zero(d1, "mu"); p2 <- p_at_zero(d2, "mu")
    },
    "fixed_opt_avg" = {
      d0 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 0)
      d1 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 1)
      d2 <- get_fixed_opt_avg(con, bf_crit, r_val, decision = 2)
      p0 <- p_at_zero(d0, "delta"); p1 <- p_at_zero(d1, "delta"); p2 <- p_at_zero(d2, "delta")
    },
    "weighted_sum" = {
      d0 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 0)
      d1 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 1)
      d2 <- get_fixed_weighted_sum(con, bf_crit, r_val, decision = 2)
      p0 <- p_at_zero(d0, "delta"); p1 <- p_at_zero(d1, "delta"); p2 <- p_at_zero(d2, "delta")
    },
    stop("Invalid df_type")
  )

  c00 * p0 + c10 * p1 + ci0 * p2
}

# NEW: theoretical Cauchy prior mass captured within [lower, upper]
cauchy_interval_mass <- function(r_val, lower = -1, upper = 1) {
  pcauchy(upper, location = 0, scale = r_val) - pcauchy(lower, location = 0, scale = r_val)
}


###########################################
############ WEIGHT COSTS #################
###########################################

extend_grid_interp <- function(x, y, r_val, target_delta,
                                grid_extent_mode = c("interval", "threshold"),
                                mass_interval = c(-1, 1),
                                tail_threshold = 1e-7) {
  grid_extent_mode <- match.arg(grid_extent_mode)
  max_x <- max(x)
  n <- length(x)
  slope <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])  # NEW: moved up, reused in loop + final extrapolation

  interval_bound <- switch(grid_extent_mode,
"interval"  = max(abs(mass_interval)),
"threshold" = FALSE
  )

# interpolate for cauchy weight decay
if (!interval_bound){
    grid <- seq(0, max_x, by = target_delta)
    step <- max_x
repeat {
      step <- step + target_delta
      w <- target_delta / (pi * r_val * (1 + (step / r_val)^2))
      cost_at_step <- y[n] + slope * (step - x[n])   
if (w * cost_at_step < tail_threshold) break # we need to account for weight
      grid <- c(grid, step)
    }
  }else{
    upper_bound <- max(max_x, interval_bound)
    grid <- seq(0, upper_bound, by = target_delta)
  }

  y_grid <- approx(x, y, xout = grid, rule = 2)$y

  # linear extrapolation based on the slope of the last two points
  extended_idx <- grid > max_x
  y_grid[extended_idx] <- y[n] + slope * (grid[extended_idx] - x[n])

  data.frame(delta = grid, cost = y_grid)
}

plot_prior_distr_cost <- function(con, bf_crit, r_val, cost_matrix,
                                   df_type = "opt_stop", target_delta = 0.05,
                                   grid_extent_mode = c("interval", "threshold"),
                                   mass_interval = c(-1, 1),
                                   tail_threshold = 1e-7,
                                   scenario_name = "default") {
  grid_extent_mode <- match.arg(grid_extent_mode)

  col <- switch(df_type,
    "opt_stop"      = my_colors[2],
    "fixed_max"     = my_colors[3],
    "fixed_opt_avg" = my_colors[4],
    "weighted_sum"  = my_colors[5],
    stop("Invalid df_type")
  )
  title_text <- switch(df_type,
    "opt_stop" = "Optional Stopping", "fixed_max" = "MAX",
    "fixed_opt_avg" = "STOP AVG", "weighted_sum" = "SAME DIST"
  )

  df <- get_cost_curve(con, bf_crit, r_val, cost_matrix, df_type)
  delta_or_mu <- if ("delta" %in% names(df)) "delta" else "mu"
  names(df)[names(df) == delta_or_mu] <- "delta"

  df <- aggregate(cost ~ delta, data = df, FUN = mean)
  df <- df[order(df$delta), ]

  df_grid <- extend_grid_interp(df$delta, df$cost, r_val, target_delta,
                                 grid_extent_mode, mass_interval, tail_threshold)

  deltas  <- df_grid$delta[df_grid$delta > 0]
  delta_0 <- df_grid$delta[df_grid$delta == 0]

  weights   <- target_delta / (pi * r_val * (1 + (deltas / r_val)^2))
  weights_0 <- target_delta / (pi * r_val * (1 + (delta_0 / r_val)^2))

  weighted_cost   <- df_grid$cost[df_grid$delta > 0]  * weights
  weighted_cost_0 <- df_grid$cost[df_grid$delta == 0] * weights_0

  deltas_plot <- c(-rev(deltas), delta_0, deltas)
  cost_plot   <- c(rev(weighted_cost), weighted_cost_0, weighted_cost)

  # diagnostics
  grid_extent   <- max(deltas)
  captured_mass <- pcauchy(grid_extent, location = 0, scale = r_val) -
                    pcauchy(-grid_extent, location = 0, scale = r_val)
  h0_point_cost <- get_h0_point_cost(con, bf_crit, r_val, cost_matrix, df_type)
  h1_total_cost <- sum(cost_plot)
  expected_cost <- 0.5 * h0_point_cost + 0.5 * h1_total_cost   # NEW: 50/50 blend across ground truths

  extent_label <- if (grid_extent_mode == "interval") {
    sprintf("interval [%.3g,%.3g]", mass_interval[1], mass_interval[2])
  } else {
    sprintf("tail threshold %.1g", tail_threshold)
  }

  print(data.frame(delta = deltas_plot, weighted_cost = cost_plot))
  cat(sprintf(
    "[%s, BF_crit=%s, r=%.3f, mode=%s] grid extent=[%.3g,%.3g] | mass captured=%.6f | cost|H0=%.4f | cost|H1=%.4f | E[cost]=%.4f\n",
    df_type, bf_crit, r_val, extent_label, -grid_extent, grid_extent, captured_mass,
    h0_point_cost, h1_total_cost, expected_cost
  ))

  # NEW: mode goes into the filename so interval/threshold don't overwrite each other
  pdf(sprintf("figures/prior-weighted-cost_%s-%s-bf%s-%s.pdf",
              df_type, scenario_name, bf_crit, grid_extent_mode))
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
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = h0_point_cost,
          col = "black", lwd = 3, lty = 2)
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
  invisible(list(
    data = data.frame(delta = deltas_plot, weighted_cost = cost_plot),
    h0_point_cost = h0_point_cost,
    h1_total_cost = h1_total_cost,
    expected_cost = expected_cost,
    interval_mass = captured_mass,
    grid_extent = grid_extent
  ))
}

# NEW: map each BF_crit to its purpose-built cost matrix
bf_crit_cost_matrices <- list(
  `3`  = cost_matrix4,
  `6`  = cost_matrix5,
  `10` = cost_matrix6
)

con <- dbConnect(duckdb(), db_file)
for (bfc in c(3, 6, 10)) {
  cm <- bf_crit_cost_matrices[[as.character(bfc)]]
  for (dtype in c("opt_stop", "fixed_max", "fixed_opt_avg", "weighted_sum")) {
    for (mode in c("interval", "threshold")) {   # NEW: both extent modes, every call
      plot_prior_distr_cost(con, bf_crit = bfc, r_val = r_vals[2],
                             cost_matrix = cm, df_type = dtype,
                             grid_extent_mode = mode,
                             mass_interval = c(-1, 1),
                             tail_threshold = 10e-5,
                             scenario_name = sprintf("sym-ind-bfcrit%s", bfc))
    }
  }
}
dbDisconnect(con)