# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("scripts/plots/define_colors.R")
# Initialize variables
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

# Convert CSV to SQLite DB
db_file <- "data/hacking-bayes.duckdb"
table_name_opt_stop <- "realistic_sym_simulation"
table_name_cauchy_fun <- "cauchy_prior"


####################
##### QUERIES ######
####################
get_opt_stop <- function(con, bf_crit, r_val, decision = 0) {
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
  fixed_max <- dbGetQuery(
    con, "
    WITH prob_table AS (
      SELECT MU, trial_count, COUNT (CASE WHEN decision = ? THEN 1 END) * 1.0 / COUNT(*) AS prob
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
      WHERE ABS(p2.mu - p.mu) < 1e-6
    )
    ORDER BY mu, trial_count",
    list(decision, r_val, bf_crit)
  )
  fixed_max
}

get_prob_by_trial <- function(con, bf_crit, r_val, decision = 1, fixed_max_df) {
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
        END AS prob_interp
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
    SELECT delta, avg_stop, prob_opt, prob_interp
    FROM interp
    ORDER BY delta",
    params = list(decision, r_val, bf_crit, decision, r_val, bf_crit)
  )
  fixed_opt_avg
}

get_fixed_weighted_sum <- function(con, bf_crit, r_val, decision = 0) {
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
  SUM(weight_component) AS weighted_sum
FROM weighted
GROUP BY mu
ORDER BY delta;
", list(r_val, bf_crit, decision, r_val, bf_crit))
  weighted_sum
}

#######################
####### PLOTTING ######
#######################

cauchy_plot <- function() {
  # read the simulation results from the file
  con <- dbConnect(duckdb(), db_file)
  # get mus
  mus <- unlist(dbGetQuery(con, "SELECT DISTINCT mu FROM cauchy_prior"))
  # average BF over repetitions
  results_df <- dbGetQuery(con, "SELECT mu,n,AVG(bf) FROM cauchy_prior GROUP BY n, mu ORDER BY mu, n")
  dbDisconnect(con)
  pdf("figures/cauchy-prior-simulation.pdf")
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
}

realistic_sim_overview_plot <- function() {
  con <- dbConnect(duckdb(), db_file)
  pdf(paste("figures/realistic-sym-decision-prob", ".pdf", sep = ""), width = 12, height = 8)
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
  base <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[1])
  )
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
  bf1 <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[1])
  )
  bf2 <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[2])
  )
  bf3 <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[3])
  )
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
  r1 <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                    AND trial_start = 2 AND trial_end = 100000
                    GROUP BY mu
                    ORDER BY mu",
    list(r_vals[1], BF_crits[1])
  )
  r2 <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                    AND trial_start = 2 AND trial_end = 100000
                    GROUP BY mu
                    ORDER BY mu",
    list(r_vals[2], BF_crits[1])
  )
  r3 <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                    AND trial_start = 2 AND trial_end = 100000
                    GROUP BY mu
                    ORDER BY mu",
    list(r_vals[3], BF_crits[1])
  )
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
  short_range <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 20 AND trial_end = 200
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[1])
  )
  long_range <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[1])
  )
  lines(long_range[[1]], long_range[[2]], col = my_colors[1], lwd = 3)
  lines(short_range[[1]], short_range[[2]], col = my_colors[2], lwd = 3)
  abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
  legend("topright",
    legend = c("2 to 100000", "20 to 200"),
    col = c(my_colors[1], my_colors[2]), lwd = 3
  )
  dev.off()
  # }
  dbDisconnect(con)
}
# Plot three BF01 curves for mus and decision probability for H0 given mu
realistic_sim_histograms <- function(bf_crit, mus, r_val) {
  con <- dbConnect(duckdb(), db_file)
  # Three Histograms for mus with H0 and H1 decisions
  custom_breaks <- seq(0, 4000, by = 1)
  r_text <- round(r_val, 3)
  for (mu in mus) {
    h0_data <- dbGetQuery(
      con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND ABS(? - mu) < 1e-6
                    AND decision = 0",
      list(r_val, bf_crit, mu)
    )
    h1_data <- dbGetQuery(
      con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND ABS(? - mu) < 1e-6
                    AND decision = 1",
      list(r_val, bf_crit, mu)
    )
    indecisive_data <- dbGetQuery(
      con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND ABS(? - mu) < 1e-6
                    AND decision = 2",
      list(r_val, bf_crit, mu)
    )
    pdf(paste("figures/realistic-sim-bf-crit-", bf_crit, "-r-", r_text, "-mu-", mu, ".pdf", sep = ""))
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
  dbDisconnect(con)
}

realistic_sim_decision_prob_curve <- function(bf_crit, r_val) {
  con <- dbConnect(duckdb(), db_file)
  # Decision probability for H0 given mu
  r_text <- round(r_val, 3)
  pdf(paste("figures/realistic-sim-decision-prob-bf-crit2-", bf_crit, "-r-", r_text, ".pdf", sep = ""), width = 14, height = 8)
  par(mar = c(5, 6, 5, 5))
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability for " * H[0] * " with " * BF[crit] * " = " * .(bf_crit) * " and " * r * " = " * .(r_text)),
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(delta),
    cex.axis = 2, cex.lab = 2, cex.main = 2
  )
  long_range <- dbGetQuery(
    con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
    list(r_vals[2], BF_crits[1])
  )
  # get the first index where the probability is smaller than 0.5
  lines(long_range[["mu"]], long_range[["prob"]], col = "black", lwd = 5)
  dev.off()
  dbDisconnect(con)
}
realistic_sim_decision_prob_curve(3, r_vals[2])

# WORK IN PROGRESS!
# db_file <- "data/results.duckdb"# WORK IN PROGRESS!
# db_file <- "data/results.duckdb"

fixed_big_mus <- c(seq(0, 0.09, 0.01), seq(0.1, 1, 0.05))
BF_crits <- c(3)
# Open connection
con <- dbConnect(duckdb(), db_file)

# Precompute all decision probability summaries
fixed_opt_avg_data <- get_fixed_opt_avg(con, BF_crits[1], r_vals[2])
fixed_weighted_sum_data <- get_fixed_weighted_sum(con, BF_crits[1], r_vals[2])
fixed_max_data <- get_fixed_max(con, BF_crits[1], r_vals[2])

# Disconnect after plots
dbDisconnect(con)

# Plot decision probability for 'H0' on y axis and sample size on the x axis
realistic_sim_fixed_size_plot <- function(fixed_opt_avg_data, fixed_weighted_sum_data, fixed_max_data) {
  # increase the size, including font size
  con <- dbConnect(duckdb(), db_file)
  for (mu in c(0.1)) {
    pdf(paste("figures/fixed-sim/realistic-fixed-decision-prob-bfcrit-", BF_crits[1], "-r-", r_vals[2], "-mu-", mu, "-comp.pdf", sep = ""), width = 12, height = 8)
    par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
    plot(0, 0,
      xlim = c(0, 1000), ylim = c(0, 1), type = "n",
      main = "",
      ylab = bquote("Decision Probability for " * H[0]), xlab = "n",
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

      # Prepare histogram data
      stop_counts_expanded <- rep(opt_stop$stop_count, opt_stop$count)
      # get maximum stop count
      max_stop_count <- max(opt_stop$stop_count)
      # Create histogram data with breaks
      hist_data <- hist(stop_counts_expanded, breaks = c(seq(0, max_stop_count, 1)), plot = FALSE)

      # Normalize histogram counts to fit within y = (0, 1)
      hist_heights <- hist_data$counts / max(hist_data$counts)

      # Draw histogram on the plot using scaled heights
      # for (i in 1:length(hist_data$counts)) {
      #   rect(
      #     xleft = hist_data$breaks[i],
      #     xright = hist_data$breaks[i + 1],
      #     ybottom = 0,
      #     ytop = hist_heights[i],
      #     col = my_colors[7], # adjustcolor(my_colors[4], alpha.f = 0.4),
      #     border = NA
      #   )
      # }
      abline(h = opt_prob_mean, col = my_colors[2], lwd = 5, lty = 2)
      # abline(v = opt_stop_mean, col = my_colors[4], lwd = 5)
      lines(fixed_size[[1]], fixed_size[[2]], col = my_colors[1], lwd = 5)
      # abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
      interp_val <- fixed_opt_avg_data[abs(fixed_opt_avg_data$delta - mu) < 1e-6, "prob_interp"]
      weighted_val <- fixed_weighted_sum_data[abs(fixed_weighted_sum_data$delta - mu) < 1e-6, "weighted_sum"]
      max_val <- fixed_max_data[abs(fixed_max_data$mu - mu) < 1e-6, "prob"]
      if (length(max_val) == 1) {
        abline(h = max_val, col = my_colors[3], lwd = 5, lty = 2)
      }
      # if (length(interp_val) == 1) {
      #   abline(h = interp_val, col = my_colors[4], lwd = 5, lty = 2)
      # }
      # if (length(weighted_val) == 1) {
      #   abline(h = weighted_val, col = my_colors[5], lwd = 5, lty = 2)
      # }
      legend("bottomright",
        legend = c("fixed n", "optional stopping", "MAX"), # "STOP AVG", "SAME DIST"),
        col = c(my_colors[1], my_colors[2], my_colors[3]), # my_colors[4], my_colors[5]),
        lwd = 5, lty = c(1, 2, 2), # 2, 2)
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
  par(mfrow = c(1, 1), mar = c(15, 1, 1, 1))
  pdf(paste("figures/realistic-fixed-sim-overview-for-bfcrit-", bf_crit, "-r-", round(r_val, 3), ".pdf", sep = ""))
  par(mar = c(5, 6, 5, 5))
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability with " * BF[crit] == .(bf_crit) * ", r =" * .(round(r_val, 3))),
    ylab = bquote("Decision Probability for " * H[0]),
    xlab = bquote("Effect size " * delta),
    cex.main = 2, cex.lab = 2, cex.axis = 2
  )
  lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd = 5)
  lines(fixed_max[["mu"]], fixed_max[["prob"]], col = my_colors[3], lwd = 5)
  lines(fixed_opt_avg[["delta"]], fixed_opt_avg[["prob_interp"]], col = my_colors[4], lwd = 5)
  lines(weighted_results[["delta"]], weighted_results[["weighted_sum"]], col = my_colors[5], lwd = 5)

  legend("topright",
    legend = c("optional stopping", "MAX", "STOP AVG", "SAME DIST"),
    col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
    lwd = 5
  )
  dev.off()
}

realistic_sim_combined_plot(3, r_vals[3])
realistic_sim_combined_plot(3, r_vals[2])
realistic_sim_combined_plot(3, r_vals[1])
realistic_sim_combined_plot(6, r_vals[2])
realistic_sim_combined_plot(10, r_vals[2])

stop_count_to_effect_size <- function(con, bf_crit, r_val) {
  # Query with rounding and correct GROUP BY
  stop_count_data <- dbGetQuery(
    con, "
      SELECT mu, AVG(stop_count) AS avg_stop_count
      FROM cauchy_sym
      WHERE bf_crit = ? AND ABS(r - ?) < 1e-6
      AND trial_start = 2 AND trial_end = 100000
      GROUP BY mu
      ORDER BY mu
    ",
    list(bf_crit, r_val)
  )
  stop_count_data
}

# plot for all rs and for all BF_crits
plot_stop_count_to_effect_size_bf <- function(bf_crits, r_val, db_file = "data/hacking-bayes.duckdb") {
  con <- dbConnect(duckdb(), db_file)
  pdf(paste("figures/stop-count-to-effect-size-r-", round(r_val, 3), ".pdf", sep = ""), width = 12, height = 8)
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 400), type = "n",
    main = bquote("Stop Count to Effect Size for " * r * " = " * .(r_val)),
    ylab = bquote("Stop Count"),
    xlab = bquote("Effect Size " * delta)
  )
  for (i in seq_along(bf_crits)) {
    stop_count_data <- stop_count_to_effect_size(con, bf_crits[i], r_val)
    # Plot the stop count data
    lines(stop_count_data$mu, stop_count_data$avg_stop_count, col = my_colors[i], lwd = 2)
  }
  dbDisconnect(con)
  # Add a legend entry for each BF_crit
  legend("topright",
    legend = c(bquote(BF[crit] == .(bf_crits[1])), bquote(BF[crit] == .(bf_crits[2])), bquote(BF[crit] == .(bf_crits[3]))),
    col = c(my_colors[1], my_colors[2], my_colors[3]), lwd = 2, cex = 0.8
  )
  dev.off()
}

# plot for all rs
plot_stop_count_to_effect_size_r <- function(bf_crit, r_vals, db_file = "data/hacking-bayes.duckdb") {
  con <- dbConnect(duckdb(), db_file)
  pdf(paste("figures/stop-count-to-effect-size-bf-", bf_crit, ".pdf", sep = ""), width = 12, height = 8)
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 50), type = "n",
    main = bquote("Stop Count to Effect Size for " * BF[crit] * " = " * .(bf_crit)),
    ylab = bquote("Stop Count"),
    xlab = bquote("Effect Size " * delta)
  )
  for (i in seq_along(r_vals)) {
    stop_count_data <- stop_count_to_effect_size(con, bf_crit, r_vals[i])
    # Plot the stop count data
    lines(stop_count_data$mu, stop_count_data$avg_stop_count, col = my_colors[i], lwd = 2)
  }
  dbDisconnect(con)
  # Add a legend entry for each r_val
  legend("topright",
    legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
    col = c(my_colors[1], my_colors[2], my_colors[3]), lwd = 2, cex = 0.8
  )
  dev.off()
}

plot_stop_count_to_effect_size_bf(c(3, 6, 10), 1 / sqrt(2))
plot_stop_count_to_effect_size_r(3, c(0.5 / sqrt(2), 1 / sqrt(2), 2 / sqrt(2)))

# plot difference between optional stopping and fixed size for the same mu
plot_fixed_vs_optional_stopping <- function(bf_crit, r_val, db_file = "data/hacking-bayes.duckdb") {
  con <- dbConnect(duckdb(), db_file)
  pdf(paste("figures/fixed-vs-optional-stopping-bf-", bf_crit, "-r-", round(r_val, 3), ".pdf", sep = ""), width = 12, height = 8)
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Fixed Size vs Optional Stopping for " * BF[crit] * " = " * .(bf_crit) * ", r = " * .(r_val)),
    ylab = bquote("Decision Probability for " * H[0]),
    xlab = bquote("Effect Size " * delta)
  )

  fixed_max <- get_fixed_max(con, bf_crit, 2 * r_val)
  fixed_opt_avg <- get_fixed_opt_avg(con, bf_crit, 2 * r_val)
  weighted_sum <- get_fixed_weighted_sum(con, bf_crit, 2 * r_val)
  opt_stop <- get_opt_stop(con, bf_crit, r_val)

  lines(fixed_max$mu, fixed_max$prob, col = my_colors[3], lwd = 2, lty = 2)
  lines(opt_stop$mu, opt_stop$prob, col = my_colors[2], lwd = 2)
  lines(fixed_opt_avg$delta, fixed_opt_avg$prob_interp, col = my_colors[4], lwd = 2, lty = 4)
  lines(weighted_sum$delta, weighted_sum$weighted_sum, col = my_colors[5], lwd = 2, lty = 5)

  legend("topright",
    legend = c("optional stopping prob", bquote(max[fixed] * " " * n), bquote(bar(n)[os] * "as fixed size"), "same distribution as fixed"),
    col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
    lwd = 2, lty = c(1, 2, 4, 5)
  )

  dev.off()

  dbDisconnect(con)
}
plot_fixed_vs_optional_stopping(3, r_vals[2])

plot_decision_sim_fixed_size <- function(bf_crit, r_val) {
  pdf_filename <- sprintf("figures/realistic-sim-fixed-size-all-decisions-bf-crit-%s-r-%.3f.pdf", bf_crit, r_val)
  pdf(pdf_filename, width = 16, height = 8)
  par(mfrow = c(1, 3), mar = c(1, 1, 1, 1))

  con <- dbConnect(duckdb(), db_file)
  fixed_max_d0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)

  dbDisconnect(con)

  for (decision in 0:2) {
    decision_text <- switch(as.character(decision),
                            "0" = "H0",
                            "1" = "H1",
                            "2" = "Indecisive")

    con <- dbConnect(duckdb(), db_file)

    # Verwende nur get_prob_by_trial() für decision != 0
    if (decision == 0) {
      fixed_max <- fixed_max_d0
    } else {
      fixed_max <- get_prob_by_trial(con, bf_crit, r_val, decision, fixed_max_df = fixed_max_d0)
    }
    print(fixed_max)
    opt_stop <- get_opt_stop(con, bf_crit, r_val, decision)
    fixed_opt_avg <- get_fixed_opt_avg(con, bf_crit, r_val, decision)
    weighted_results <- get_fixed_weighted_sum(con, bf_crit, r_val, decision)
    dbDisconnect(con)

    par(mar = c(5, 6, 5, 5))
    plot(0, 0,
         xlim = c(0, 1), ylim = c(0, 1), type = "n",
         main = decision_text,
         ylab = bquote("Decision Probability for " * .(decision_text)),
         xlab = bquote("Effect size " * delta),
         cex.main = 2, cex.lab = 2, cex.axis = 2)

    if (nrow(opt_stop) > 0) lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd = 5)
    if (nrow(fixed_max) > 0) lines(fixed_max[["mu"]], fixed_max[["prob"]], col = my_colors[3], lwd = 5)
    if (nrow(fixed_opt_avg) > 0) lines(fixed_opt_avg[["delta"]], fixed_opt_avg[["prob_interp"]], col = my_colors[4], lwd = 5)
    if (nrow(weighted_results) > 0) lines(weighted_results[["delta"]], weighted_results[["weighted_sum"]], col = my_colors[5], lwd = 5)

    legend("topright",
           legend = c("optional stopping", "MAX", "STOP AVG", "SAME DIST"),
           col = c(my_colors[2], my_colors[3], my_colors[4], my_colors[5]),
           lwd = 5)
  }

  dev.off()
}

plot_decision_sim_fixed_size(3, r_vals[2])