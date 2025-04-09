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
                    AND mu = ?
                    AND decision = 0",
            list(r_val, bf_crit, mu)
        )
        h1_data <- dbGetQuery(
            con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND mu = ?
                    AND decision = 1",
            list(r_val, bf_crit, mu)
        )
        indecisive_data <- dbGetQuery(
            con, "SELECT stop_count
                    FROM cauchy_sym
                    WHERE ABS(? - r) < 1e-6
                    AND bf_crit = ?
                    AND mu = ?
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
    pdf(paste("figures/realistic-sim-decision-prob-bf-crit2-", bf_crit, "-r-", r_text, ".pdf", sep = ""))
    plot(0, 0,
        xlim = c(0, 1), ylim = c(0, 1), type = "n",
        main = bquote("Decision probability for symmetrical rule with " * BF[crit] * " = " * .(bf_crit) * " and " * r * " = " * .(r_text)),
        ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
    )
    r_prob <- dbGetQuery(
        con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                    FROM cauchy_sym
                    WHERE bf_crit = ? AND ABS(? - r) < 1e-6
                    GROUP BY mu
                    ORDER BY mu",
        list(bf_crit, r_val)
    )
    is_ph0_smaller_50 <- r_prob < 0.5
    # get the first index where the probability is smaller than 0.5
    first_index <- which(is_ph0_smaller_50)[1]
    print(paste("P(H0) > 0.5 for mu =", big_sim_mus[first_index]))
    lines(r_prob[["mu"]], r_prob[["prob"]], col = "black", lwd = 2)
    dev.off()
    dbDisconnect(con)
}

# WORK IN PROGRESS!
# db_file <- "data/results.duckdb"
fixed_big_mus <- seq(0, 1, 0.05)
# Plot decision probability for 'H0' on y axis and sample size on the x axis
realistic_sim_fixed_size_plot <- function() {
    con <- dbConnect(duckdb(), db_file)
    for (mu in fixed_big_mus) {
        pdf(paste("figures/fixed-sim/realistic-fixed-size-decision-prob-", mu, ".pdf", sep = ""), width = 12, height = 8)
        # base case
        plot(0, 0,
            xlim = c(0, 1000), ylim = c(0, 1), type = "n",
            main = "",
            ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(delta),
            cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
        )
        for (BF_crit in BF_crits) {
            fixed_size <- dbGetQuery(
                con, "SELECT trial_count, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob
                FROM cauchy_sym_fixed_size
                WHERE mu = ? AND ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 1000
                GROUP BY trial_count
                ORDER BY trial_count",
                list(mu, r_vals[2], BF_crit)
            )
            opt_stop <- dbGetQuery(
                con, "SELECT stop_count, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob,
                COUNT(stop_count) AS count
                FROM cauchy_sym
                WHERE mu = ? AND ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY stop_count
                ORDER BY stop_count",
                list(mu, r_vals[2], BF_crit)
            )
            opt_stop_mean <- sum(opt_stop[[1]] * opt_stop[[3]] / 20000)
            opopt_stop_sd <- sd(opt_stop[[1]])
            opt_prob_mean <- sum(opt_stop[[2]] * opt_stop[[3]] / 20000)
            opt_prob_sd <- sd(opt_stop[[2]])

            points(opt_stop_mean, opt_prob_mean, pch = 4)
            abline(h = opt_prob_mean, col = my_colors[1])
            lines(fixed_size[[1]], fixed_size[[2]], col = my_colors[1], lwd = 3)
            abline(h = 1 / BF_crits[1], col = my_colors[1], lty = 2, lwd = 3)
        }
        dev.off()
    }
    dbDisconnect(con)
}
realistic_sim_fixed_size_plot()

# Plot maximum probability of 'H0' for each mu

realistic_sim_fixed_max <- function() {
    con <- dbConnect(duckdb(), db_file)
    # Decision probability for H0 given mu
    pdf(paste("figures/realistic-sim-fixed-max-bf10.pdf", sep = ""))
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
                    AND trial_start = 2 AND trial_end = 1000
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
        list(r_vals[2], BF_crits[3])
    )
    opt_stop <- dbGetQuery(
        con, "SELECT mu, COUNT (CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS prob,
                COUNT(stop_count) AS count
                FROM cauchy_sym
                WHERE ABS(? - r) < 1e-6 AND bf_crit = ?
                AND trial_start = 2 AND trial_end = 100000
                GROUP BY mu
                ORDER BY mu",
        list(r_vals[2], BF_crits[3])
    )
    # get the first index where the probability is smaller than 0.5
    lines(fixed_max[["mu"]], fixed_max[["prob"]], col = "black", lwd = 2)
    lines(opt_stop[["mu"]], opt_stop[["prob"]], col = my_colors[2], lwd =2)
    dev.off()
    dbDisconnect(con)
}
realistic_sim_fixed_max()
