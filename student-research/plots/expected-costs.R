# Set working directory
setwd(".")

# Clear workspace
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
db_large_effects_file <- "data/large_effects.duckdb"


###########################################
############### QUERIES  ##################
###########################################

get_opt_stop <- function(con_large, bf_crit, r_val, decision = 0) {
    dbExecute(con_large, "SET max_expression_depth TO 10000")
    dbGetQuery(
        con_large, "
    SELECT mu, COUNT (CASE WHEN decision = ? THEN 1.0 END) * 1.0 / COUNT(*) AS prob,
           AVG(stop_count) AS mean_count
    FROM cauchy_sym
    WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
      AND trial_start = 2 AND trial_end = 100000
    GROUP BY mu
    ORDER BY mu",
        list(decision, r_val, bf_crit)
    )
}

get_fixed_max <- function(con_large, bf_crit, r_val, decision = 0) {
    dbExecute(con_large, "SET max_expression_depth TO 10000")
    dbGetQuery(
        con_large, "
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
}

get_prob_by_trial <- function(con_large, bf_crit, r_val, decision = 1, fixed_max_df) {
    dbExecute(con_large, "SET max_expression_depth TO 10000")
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

    dbGetQuery(con_large, sql, list(decision, r_val, bf_crit))
}
get_fixed_stop_avg <- function(con, bf_crit, r_val, decision = 0) {
    dbExecute(con, "SET max_expression_depth TO 10000")
    fixed_stop_avg <- dbGetQuery(con, "
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
    fixed_stop_avg
}

get_fixed_same_dist <- function(con, bf_crit, r_val, decision = 0) {
    dbExecute(con, "SET max_expression_depth TO 10000")
    same_dist <- dbGetQuery(con, "
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
    same_dist
}



#####################################
####### Expected cost summaries #########
#####################################

# cost_matrix is a 2x3 list
# ====================================================================
# REFACTORED: Summary function for expected costs (original plot)
# ====================================================================
summary_expected_costs_original <- function(
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
    db_file = "data/hacking-bayes.duckdb") {
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
    fa0 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 0)
    fa1 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 1)
    fa2 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 2)
    fa_delta <- fa0$delta
    fa_cost <- ec(
        fa0$prob,
        reagrid(fa1$delta, fa1$prob, fa_delta),
        reagrid(fa2$delta, fa2$prob, fa_delta),
        fa_delta
    )

    # Fixed SAME DIST — uses delta as the x-axis
    ws0 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 0)
    ws1 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 1)
    ws2 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 2)
    ws_delta <- ws0$delta
    ws_cost <- ec(
        ws0$prob,
        reagrid(ws1$delta, ws1$prob, ws_delta),
        reagrid(ws2$delta, ws2$prob, ws_delta),
        ws_delta
    )

    dbDisconnect(con)

    # Return summary data structure
    list(
        bf_crit = bf_crit,
        r_val = r_val,
        cost_matrix = cost_matrix,
        c00 = c00, c10 = c10, ci0 = ci0,
        c01 = c01, c11 = c11, ci1 = ci1,
        os_mu = os_mu, os_cost = os_cost,
        fm_mu = fm_mu, fm_cost = fm_cost,
        fa_delta = fa_delta, fa_cost = fa_cost,
        ws_delta = ws_delta, ws_cost = ws_cost
    )
    }


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
        0, 1, 1 / 7,
        1, 0, 1 / 7
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)

cost_matrix6 <- matrix(
    c(
        0, 1, 1 / 11,
        1, 0, 1 / 11
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
        is_zero <- abs(delta_vals) < 1e-9 # FIXED: was `is_h0`, undefined at point of use
        ifelse(is_zero,
            c01 * p0 + ci1 * p2, # delta = 0: H1 ground truth; only dec-H0 and indecisive costs
            c01 * p0 + c11 * p1 + ci1 * p2 # delta != 0: unchanged, full H1-true blend
        )
    }
    reagrid <- function(src_x, src_y, target_x) {
        approx(src_x, src_y, xout = target_x, rule = 2)$y
    }
    # delta_or_mu <- if ("delta" %in% names(df)) "delta" else "mu"
    # names(df)[names(df) == delta_or_mu] <- "delta"

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
        "stop_avg" = {
            d0 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 0)
            d1 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 1)
            d2 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 2)
            x <- d0$delta
            data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$delta, d1$prob, x), reagrid(d2$delta, d2$prob, x), x))
        },
        "same_dist" = {
            d0 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 0)
            d1 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 1)
            d2 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 2)
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
        if (length(idx) == 0) {
            return(0)
        }
        df$prob[idx[1]]
    }

    switch(df_type,
        "opt_stop" = {
            d0 <- get_opt_stop(con, bf_crit, r_val, decision = 0)
            d1 <- get_opt_stop(con, bf_crit, r_val, decision = 1)
            d2 <- get_opt_stop(con, bf_crit, r_val, decision = 2)
            p0 <- p_at_zero(d0, "mu")
            p1 <- p_at_zero(d1, "mu")
            p2 <- p_at_zero(d2, "mu")
        },
        "fixed_max" = {
            d0 <- get_fixed_max(con, bf_crit, r_val, decision = 0)
            d1 <- get_prob_by_trial(con, bf_crit, r_val, decision = 1, fixed_max_df = d0)
            d2 <- get_prob_by_trial(con, bf_crit, r_val, decision = 2, fixed_max_df = d0)
            p0 <- p_at_zero(d0, "mu")
            p1 <- p_at_zero(d1, "mu")
            p2 <- p_at_zero(d2, "mu")
        },
        "stop_avg" = {
            d0 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 0)
            d1 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 1)
            d2 <- get_fixed_stop_avg(con, bf_crit, r_val, decision = 2)
            p0 <- p_at_zero(d0, "delta")
            p1 <- p_at_zero(d1, "delta")
            p2 <- p_at_zero(d2, "delta")
        },
        "same_dist" = {
            d0 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 0)
            d1 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 1)
            d2 <- get_fixed_same_dist(con, bf_crit, r_val, decision = 2)
            p0 <- p_at_zero(d0, "delta")
            p1 <- p_at_zero(d1, "delta")
            p2 <- p_at_zero(d2, "delta")
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
    slope <- (y[n] - y[n - 1]) / (x[n] - x[n - 1]) # NEW: moved up, reused in loop + final extrapolation

    interval_bound <- switch(grid_extent_mode,
        "interval"  = max(abs(mass_interval)),
        "threshold" = FALSE
    )

    # interpolate for cauchy weight decay
    if (!interval_bound) {
        grid <- seq(0, max_x, by = target_delta)
        step <- max_x
        repeat {
            step <- step + target_delta
            w <- target_delta / (pi * r_val * (1 + (step / r_val)^2))
            cost_at_step <- y[n] + slope * (step - x[n])
            if (w * cost_at_step < tail_threshold) break # we need to account for weight
            grid <- c(grid, step)
        }
    } else {
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

    df <- get_cost_curve(con, bf_crit, r_val, cost_matrix, df_type)
    delta_or_mu <- if ("delta" %in% names(df)) "delta" else "mu"
    names(df)[names(df) == delta_or_mu] <- "delta"

    df <- aggregate(cost ~ delta, data = df, FUN = mean)
    df <- df[order(df$delta), ]

    df_grid <- extend_grid_interp(
        df$delta, df$cost, r_val, target_delta,
        grid_extent_mode, mass_interval, tail_threshold
    )

    deltas <- df_grid$delta[df_grid$delta > 0]
    delta_0 <- df_grid$delta[df_grid$delta == 0]

    # print(data.frame(df_grid$delta, df_grid$cost))

    weights <- target_delta / (pi * r_val * (1 + (deltas / r_val)^2))
    weights_0 <- target_delta / (pi * r_val * (1 + (delta_0 / r_val)^2))

    cauchy_mass <- 2 * sum(weights) + sum(weights_0)

    weighted_cost <- df_grid$cost[df_grid$delta > 0] * weights
    weighted_cost_0 <- df_grid$cost[df_grid$delta == 0] * weights_0

    deltas_plot <- c(-rev(deltas), delta_0, deltas)
    cost_plot <- c(rev(weighted_cost), weighted_cost_0, weighted_cost)

    # diagnostics
    grid_extent <- max(deltas)
    captured_mass <- pcauchy(grid_extent, location = 0, scale = r_val) -
        pcauchy(-grid_extent, location = 0, scale = r_val)
    h0_point_cost <- get_h0_point_cost(con, bf_crit, r_val, cost_matrix, df_type)
    h1_total_cost <- sum(cost_plot)
    expected_cost <- 0.5 * h0_point_cost + 0.5 * h1_total_cost # NEW: 50/50 blend across ground truths

    extent_label <- if (grid_extent_mode == "interval") {
        sprintf("interval [%.3g,%.3g]", mass_interval[1], mass_interval[2])
    } else {
        sprintf("tail threshold %.1g", tail_threshold)
    }

    # print(data.frame(delta = deltas_plot, weighted_cost = cost_plot))
    cat(sprintf(
        "[%s, BF_crit=%s, r=%.3f, mode=%s] grid extent=[%.3g,%.3g] | mass captured=%.6f |weight_mass=%.6f| cost|H0=%.4f | cost|H1=%.4f | E[cost]=%.4f\n",
        df_type, bf_crit, r_val, extent_label, -grid_extent, grid_extent, captured_mass, cauchy_mass,
        h0_point_cost, h1_total_cost, expected_cost
    ))

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
    invisible(list(
        data = data.frame(delta = deltas_plot, weighted_cost = cost_plot),
        h0_point_cost = h0_point_cost,
        h1_total_cost = h1_total_cost,
        expected_cost = expected_cost,
        interval_mass = captured_mass
    ))
}

# NEW: map each BF_crit to its purpose-built cost matrix
bf_crit_cost_matrices <- list(
    `3`  = cost_matrix4,
    `6`  = cost_matrix5,
    `10` = cost_matrix6
)

con <- dbConnect(duckdb(), db_file)
for (bfc in c(3)) { # , 6, 10)) {
    cm <- bf_crit_cost_matrices[[as.character(bfc)]]
    for (dtype in c("opt_stop", "stop_avg", "same_dist")) {
        for (mode in c("interval", "threshold")) { # NEW: both extent modes, every call
            plot_prior_distr_cost(con,
                bf_crit = bfc, r_val = r_vals[2],
                cost_matrix = cm, df_type = dtype,
                grid_extent_mode = mode,
                mass_interval = c(-1, 1),
                tail_threshold = 10e-5,
                scenario_name = sprintf("sym-ind-bfcrit%s", bfc)
            )
        }
    }
}
dbDisconnect(con)

# stop avg
# 37          1.80   0.08142341
# 38          1.85   0.07800663
# 39          1.90   0.07458985
# 40          1.95   0.07117308
# 41          2.00   0.06775630
# 42          2.05   0.06433952
# 43          2.10   0.06092274
# 44          2.15   0.05750597
# 45          2.20   0.05408919
# 46          2.25   0.05067241

###########################################
###### WORST-CASE BINNED WEIGHT COSTS ######
###########################################
#
#
#   Bin 1: [0,        target_delta)    -> Kosten am Punkt delta = 0
#   Bin 2: [target_delta, 2*target_delta) -> Kosten am Punkt delta = target_delta
#   Bin 3: [2*target_delta, 3*target_delta) -> Kosten am Punkt delta = 2*target_delta
#   ...
calculate_bins <- function(deltas, costs, target_delta, r_val,
                           mass_interval = c(-1, 1), bound = "upper") {
    cost_fun <- function(cost1, cost2, bound = "upper") {
        if (bound == "upper") cost1 else if (bound == "lower") cost2 else stop("invalid bound option")
    }
    ord <- order(deltas)
    sdeltas <- deltas[ord]
    scosts <- costs[ord]
    # print(data.frame(deltas = sdeltas, costs = scosts))
    # Duplikate in x zusammenfassen (Mittelwert), approx() braucht eindeutige x
    if (any(duplicated(sdeltas))) {
        stop("duplicated deltas in calculate_bins() found")
    }
    n <- length(sdeltas)
    max_x <- sdeltas[n]

    max_delta <- max(max_x, max(abs(mass_interval)))

    knots <- sort(unique(c(0, sdeltas[sdeltas > 0 & sdeltas < max_delta], max_delta)))
    edges <- sort(unique(round(c(knots), 8)))
    left_e <- edges[-length(edges)]
    right_e <- edges[-1]


    idx_left <- findInterval(left_e, sdeltas)
    idx_right <- findInterval(right_e, sdeltas)

    cost_left <- scosts[pmax(idx_left, 1)]
    cost_right <- scosts[pmax(idx_right, 1)]

    cost_worst <- mapply(cost_fun, cost_left, cost_right,
        MoreArgs = list(bound = bound)
    )
    # grid unterteilen


    bin_width <- right_e - left_e

    # Masse bzw. Density der Cauchy verteilung anhand der kosten berechnen
    # dens_left <- 1 / (pi * r_val * (1 + (left_e / r_val)^2))
    # dens_right <- 1 / (pi * r_val * (1 + (right_e / r_val)^2))
    dens <- pcauchy(right_e, location = 0, scale = r_val) - pcauchy(left_e, location = 0, scale = r_val)
    # mass_left <- bin_width * dens_left
    # mass_right <- bin_width * dens_right

    bin_mass <- mapply(cost_fun, dens, dens, # mass_left, mass_right,
        MoreArgs = list(bound = bound)
    )

    df <- data.frame(
        left = left_e, right = right_e,
        cost = cost_worst, mass = bin_mass
    )

    # print(df)
    df
}

get_large_effects_cost_points <- function(con_large, bf_crit, r_val, cost_matrix, df_type = "opt_stop") {
    c01 <- cost_matrix["H1_true", "dec_H0"]
    c11 <- cost_matrix["H1_true", "dec_H1"]
    ci1 <- cost_matrix["H1_true", "dec_ind"]

    ec <- function(p0, p1, p2, delta_vals) {
        is_zero <- abs(delta_vals) < 1e-9
        ifelse(is_zero,
            c01 * p0 + ci1 * p2,
            c01 * p0 + c11 * p1 + ci1 * p2
        )
    }
    reagrid <- function(src_x, src_y, target_x) approx(src_x, src_y, xout = target_x, rule = 2)$y

    switch(df_type,
        "opt_stop" = {
            d0 <- get_opt_stop(con_large, bf_crit, r_val, decision = 0)
            d1 <- get_opt_stop(con_large, bf_crit, r_val, decision = 1)
            d2 <- get_opt_stop(con_large, bf_crit, r_val, decision = 2)
            x <- d0$mu
            data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$mu, d1$prob, x), reagrid(d2$mu, d2$prob, x), x))
        },
        "fixed_max" = {
            d0 <- get_fixed_max(con_large, bf_crit, r_val, decision = 0)
            d1 <- get_prob_by_trial(con_large, bf_crit, r_val, decision = 1, fixed_max_df = d0)
            d2 <- get_prob_by_trial(con_large, bf_crit, r_val, decision = 2, fixed_max_df = d0)
            x <- d0$mu
            data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$mu, d1$prob, x), reagrid(d2$mu, d2$prob, x), x))
        },
        "stop_avg" = {
            d0 <- get_fixed_stop_avg(con_large, bf_crit, r_val, decision = 0)
            d1 <- get_fixed_stop_avg(con_large, bf_crit, r_val, decision = 1)
            d2 <- get_fixed_stop_avg(con_large, bf_crit, r_val, decision = 2)
            x <- d0$delta
            # print(d0)
            data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$delta, d1$prob, x), reagrid(d2$delta, d2$prob, x), x))
        },
        "same_dist" = {
            d0 <- get_fixed_same_dist(con_large, bf_crit, r_val, decision = 0)
            d1 <- get_fixed_same_dist(con_large, bf_crit, r_val, decision = 1)
            d2 <- get_fixed_same_dist(con_large, bf_crit, r_val, decision = 2)
            x <- d0$delta
            data.frame(delta = x, cost = ec(d0$prob, reagrid(d1$delta, d1$prob, x), reagrid(d2$delta, d2$prob, x), x))
        }
    )
}

    # Analog zu plot_prior_distr_cost, aber mit worst-case Bin-Gewichtung.
    summary_expected_costs <- function(con, bf_crit, r_val, bound = "upper", cost_matrix,
                                    df_type = "opt_stop", target_delta = 0.05,
                                    mass_interval = c(-1, 1),
                                    scenario_name = "default",
                                    con_large = NULL) {
        ## Check if actually neeeded

        df <- get_cost_curve(con, bf_crit, r_val, cost_matrix, df_type)
        delta_or_mu <- if ("delta" %in% names(df)) "delta" else "mu"
        names(df)[names(df) == delta_or_mu] <- "delta"
        df <- aggregate(cost ~ delta, data = df, FUN = mean)

        # NEU: grosse Effektstaerken (delta=2,3) aus large_effects.duckdb einmischen,
        # falls eine con_large uebergeben wurde. calculate_bins() erkennt
        # die entstehende Luecke (z.B. 1 -> 2) automatisch und baut dort einen
        # entsprechend grossen Bin, siehe weiter unten.
        if (!is.null(con_large)) {
            df_large <- get_large_effects_cost_points(con_large, bf_crit, r_val, cost_matrix, df_type)
            df <- rbind(df, df_large)
        }

        df <- df[order(df$delta), ]


        bins <- calculate_bins(df$delta, df$cost, target_delta, r_val, mass_interval, bound = bound)

        captured_mass <- pcauchy(mass_interval[2], location = 0, scale = r_val) -
            pcauchy(mass_interval[1], location = 0, scale = r_val)
        cauchy_mass <- 2 * sum(bins$mass)
        residual_mass <- 1 - cauchy_mass
        last_cost <- tail(df$cost, n = 1)
        if (bound == "upper") residual_cost <- residual_mass * last_cost else if (bound == "lower") residual_cost <- 0 else stop("bound is not a valid option")


        weighted_cost <- bins$cost * bins$mass
        h1_total_cost <- 2 * sum(weighted_cost) + residual_cost
        h0_point_cost <- get_h0_point_cost(con, bf_crit, r_val, cost_matrix, df_type)
        expected_cost <- 0.5 * h0_point_cost + 0.5 * h1_total_cost

        extent_label <- sprintf("interval [%.3g,%.3g]", mass_interval[1], mass_interval[2])

        cat(sprintf(
            "%s, %s case, BF_crit=%s, r=%.3f: \n | theoretical mass captured=%.6f | mass captured =%.6f | residual mass=%.6f \n | residual cost=%.6f |cost|H0=%.4f | cost|H1=%.4f | E[cost]=%.4f\n",
            df_type, bound, bf_crit, r_val, captured_mass, cauchy_mass,
            residual_mass, residual_cost,
            h0_point_cost, h1_total_cost, expected_cost
        ))

        list(
            bf_crit = bf_crit,
            r_val = r_val,
            df_type = df_type,
            bound = bound,
            bins = bins,
            h0_point_cost = h0_point_cost,
            h1_total_cost = h1_total_cost,
            expected_cost = expected_cost,
            interval_mass = captured_mass,
            residual_mass = residual_mass
        )
    }

plot_cost_bound <- function(summary) {
    #####
    ##### EXTRACT FOR REFACTORING AND PLOT GIVEN SUMMARY
    #####
    # Fix color an title text depending on dataframe type

    ####!! Sort for df_type and bound !!####

    df_type <- summary$df_type
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
    bf_crit <- summary$bf_crit
    r_val <- summary$r_val
    bound <- summary$bound
    h0_point_cost <- summary$h0_point_cost
    h1_total_cost <- summary$h1_total_cost
    expected_cost <- summary$expected_cost
    ## Check if really needed
    ## interval_mass <- summary$interval_mass
    ## residual_mass <- summary$residual_mass

    # Get bins, costs and the correct symmetrical order
    # of the bins across all effect sizes
    bins <- summary$bins
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

con <- dbConnect(duckdb(), db_file)
con_large <- dbConnect(duckdb(), db_large_effects_file)

all_summaries   <- list() # one summary for all results
per_bf_summaries <- list() # one summary per bf crit

for (bfc in c(3, 6, 10)) {
  cm <- bf_crit_cost_matrices[[as.character(bfc)]]
  bf_key <- as.character(bfc)
  per_bf_summaries[[bf_key]] <- list()

  for (dtype in c("opt_stop", "stop_avg", "same_dist")) {
    for (bound in c("upper", "lower")) {

      summary <- summary_expected_costs(
        con,
        bf_crit = bfc, r_val = r_vals[2],
        bound = bound,
        cost_matrix = cm, df_type = dtype,
        mass_interval = c(-3, 3),
        con_large = con_large
      )
      summary_df <- as.data.frame(summary)

      all_summaries[[length(all_summaries) + 1]] <- summary_df
      per_bf_summaries[[bf_key]][[length(per_bf_summaries[[bf_key]]) + 1]] <- summary_df
    }
  }
  df_bf <- do.call(rbind, per_bf_summaries[[bf_key]])
  write.table(
    df_bf,
    sprintf("student-research/summaries/expected_costs_bf-%s.txt", bfc),
    sep = "\t", row.names = FALSE, col.names = TRUE
  )
}

df_all <- do.call(rbind, all_summaries)
write.table(
  df_all,
  "student-research/summaries/expected_costs_all.txt",
  sep = "\t", row.names = FALSE, col.names = TRUE
)
dbDisconnect(con)
dbDisconnect(con_large, shutdown = TRUE)

# FIX: Read table and function for plot creation
for (bfc in c(3, 6, 10)) {
    cm <- bf_crit_cost_matrices[[as.character(bfc)]]
    #for (dtype in c("opt_stop", "stop_avg", "same_dist")) {
        #for (bound in c("upper", "lower")) {
            summary_title <- sprintf("student-research/summaries/expected_costs_bf-%s.txt", bfc)
            summary <- read.table(summary_title)
        #}
    #}
}

summary_all <- read.table(
  "student-research/summaries/expected_costs_all.txt",
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, quote = "\""
)

## overview plot for everything

plot_overview_expected_costs <- function(summary,
                                          file = "student-research/summaries/expected_costs_overview.pdf") {

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

plot_overview_expected_costs(summary_all)