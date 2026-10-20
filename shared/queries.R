# Shared DuckDB queries for the optional stopping (cauchy_sym) and the
# fixed sample size (cauchy_sym_fixed_size) simulations.
# All functions take an open connection and return a data.frame.
# Used by lab-project/summarise/fixed-size.R and student-research/summarise/expected-costs.R.

# Decision probability and mean stop count per mu (optional stopping)
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

# For every mu the fixed sample size with the highest decision probability (MAX)
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

# Decision probability at the (mu, trial_count) combinations chosen in fixed_max_df
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

# Fixed sample size decision probability interpolated at the mean stop count (STOP AVG)
get_fixed_stop_avg <- function(con, bf_crit, r_val, decision = 0) {
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

# Fixed sample size decision probability weighted with the stop count distribution (SAME DIST)
get_fixed_same_dist <- function(con, bf_crit, r_val, decision = 0) {
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
