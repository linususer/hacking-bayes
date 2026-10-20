# Summary functions for the expected cost analysis: they query the simulation
# databases (through shared/queries.R) and return the cost summaries.
# Sourced by student-research/plots/expected-costs.R.

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

# Summary: prior-weighted expected cost per effect size (one row per plotted point)
# together with the point cost given H0, the total cost given H1, the expected cost
# and the captured Cauchy mass.
prior_distr_cost_summary <- function(con, bf_crit, r_val, cost_matrix,
                                     df_type = "opt_stop", target_delta = 0.05,
                                     grid_extent_mode = c("interval", "threshold"),
                                     mass_interval = c(-1, 1),
                                     tail_threshold = 1e-7,
                                     scenario_name = "default",
                                     summary_path = NULL) {
    grid_extent_mode <- match.arg(grid_extent_mode)

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

    summary_df <- data.frame(
        bf_crit = bf_crit, r_val = r_val, df_type = df_type,
        scenario_name = scenario_name, grid_extent_mode = grid_extent_mode,
        delta = deltas_plot, weighted_cost = cost_plot,
        h0_point_cost = h0_point_cost, h1_total_cost = h1_total_cost,
        expected_cost = expected_cost, interval_mass = captured_mass
    )

    if (!is.null(summary_path)) {
        dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
        write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
        cat("=> Saved prior-weighted cost summary to:", summary_path, "\n")
    }

    summary_df
}

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
