# Set working directory
setwd("/home/linus/git/hacking-bayes")

# Clear workspace
rm(list = ls())
gc()
library(RSQLite)
library(inborutils)
library(dplyr)
source("scripts/plots/define_colors.R")
# Initialize variables
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

# Convert CSV to SQLite DB
db_file <- "data/realistic_sym_simulation.db"
csv_file <- "data/realistic_sym_simulation.csv"
table_name_opt_stop <- "realistic_sym_simulation"
table_name_cauchy_fun <- "cauchy_fun"

# # If tables with corresponding csv files do not exist you have to
# # convert them with this code snippet into the db.
# # You have to do this manually! So if in need after new calculation
# # uncomment this!

# inborutils::csv_to_sqlite(
#     csv_file = csv_file,
#     sqlite_file = db_file,
#     table_name = table_name_opt_stop,
#     pre_process_size = 1000,
#     chunk_size = 50000,
#     show_progress_bar = TRUE
# )
# inborutils::csv_to_sqlite(
#     csv_file = "data/cauchy_prior_simulation.csv",
#     sqlite_file = db_file,
#     table_name = table_name_cauchy_fun,
#     pre_process_size = 1000,
#     chunk_size = 50000,
#     show_progress_bar = TRUE
# )

cauchy_plot <- function() {
    # read the simulation results from the file
    con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    results_df <- tbl(con, table_name_cauchy_fun) %>% collect()
    # get mus
    mus <- unique(results_df$mu)
    # average BF over repetitions
    results_df <- results_df %>%
        group_by(n, mu) %>%
        summarise(bf = mean(bf))
    ## take the log
    # results_df$bf <- log(results_df$bf)
    # plot the results
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
    dbDisconnect(con)
}

realistic_sim_overview_plot <- function() {
    con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    sym <- tbl(con, table_name_opt_stop)
    # Load results # does not work because file is too big
    # sym <- read.csv("data/realistic_sym_simulation.csv")
    # Plot the decision probability for H0 given mu
    for (bf_crit in BF_crits) {
        pdf(paste("figures/realistic-sym-decision-prob-bf-crit-", bf_crit, ".pdf", sep = ""))
        plot(0, 0,
            xlim = c(0, 1), ylim = c(0, 1), type = "n",
            main = "Decision probability for symmetrical rule in realistic setting",
            ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
        )
        h0_data <- sym %>% filter(Decision_H0_H1 == 0)
        r1 <- h0_data %>%
            filter(abs(r - (0.5 / sqrt(2))) < 1e-6) %>%
            filter(BF_crit == !!bf_crit)
        r2 <- h0_data %>%
            filter(abs(r - (1 / sqrt(2))) < 1e-6) %>%
            filter(BF_crit == !!bf_crit)
        r3 <- h0_data %>%
            filter(abs(r - (2 / sqrt(2))) < 1e-6) %>%
            filter(BF_crit == !!bf_crit)
        r1_prob <- sapply(big_sim_mus, function(mu) {
            r1 %>%
                filter(Mu == mu) %>%
                count() %>%
                pull(n) / repetitions
        })
        r2_prob <- sapply(big_sim_mus, function(mu) {
            r2 %>%
                filter(Mu == mu) %>%
                count() %>%
                pull(n) / repetitions
        })
        r3_prob <- sapply(big_sim_mus, function(mu) {
            r3 %>%
                filter(Mu == mu) %>%
                count() %>%
                pull(n) / repetitions
        })
        lines(big_sim_mus, r1_prob, col = "black", lwd = 2)
        lines(big_sim_mus, r2_prob, col = "orange", lwd = 2)
        lines(big_sim_mus, r3_prob, col = "#009980", lwd = 2)
        legend("topright",
            legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
            col = c("black", "orange", "#009980"), lwd = 2
        )
        dev.off()
    }
    dbDisconnect(con)
}
# Plot three BF01 curves for mus and decision probability for H0 given mu
realistic_sim_histograms <- function(bf_crit, mus, r_val) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    sym <- tbl(con, table_name_opt_stop)
    # Three Histograms for mus with H0 and H1 decisions
    custom_breaks <- seq(0, 4000, by = 1)
    r_text <- round(r_val, 3)
    for (mu in mus) {
        h0_data <- sym %>%
            filter(Decision_H0_H1 == 0) %>%
            filter(Mu == mu) %>%
            filter(BF_crit == bf_crit) %>%
            filter(abs(r - r_val) < 1e-6)
        h1_data <- sym %>%
            filter(Decision_H0_H1 == 1) %>%
            filter(Mu == mu) %>%
            filter(BF_crit == bf_crit) %>%
            filter(abs(r - r_val) < 1e-6)
        pdf(paste("figures/realistic-sim-bf-crit-", bf_crit, "-r-", r_text, "-mu-", mu, ".pdf", sep = ""))
        par(mfrow = c(2, 1))
        par(mar = c(0, 5, 3, 3))
        # print(h0_data %>% select(Stop_Count) %>% pull())
        hist(h0_data %>% select(Stop_Count) %>% pull(),
            main = bquote("Simulation with 20000 repetitions for " * mu * " = " * .(mu) * ", " * BF[crit] * " = " * .(bf_crit) * ", " * r * " = " * .(r_text)),
            xlim = c(0, 60), ylim = c(0, 4000),
            xlab = "",
            ylab = bquote("Decision Count (" * H[0] * ")"),
            xaxt = "n", las = 1,
            col = "#59B3E6", breaks = custom_breaks
        )
        legend("topright", legend = c(bquote(H[0]), bquote(H[1])), fill = c("#59B3E6", "#CD1076"))
        par(mar = c(5, 5, 0, 3))
        hist(h1_data %>% select(Stop_Count) %>% pull(),
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
    con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    sym <- tbl(con, table_name_opt_stop)
    # Decision probability for H0 given mu
    r_text <- round(r_val, 3)
    pdf(paste("figures/realistic-sim-decision-prob-bf-crit-", bf_crit, "-r-", r_text, ".pdf", sep = ""))
    plot(0, 0,
        xlim = c(0, 1), ylim = c(0, 1), type = "n",
        main = bquote("Decision probability for symmetrical rule with " * BF[crit] * " = " * .(bf_crit) * " and " * r * " = " * .(r_text)),
        ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
    )
    h0_data <- sym %>%
        filter(Decision_H0_H1 == 0) %>%
        filter(BF_crit == bf_crit) %>%
        filter(abs(r - r_val) < 1e-6)
    r1_prob <- sapply(big_sim_mus, function(mu) {
        h0_data %>%
            filter(Mu == mu) %>%
            count() %>%
            pull(n) / repetitions
    })
    is_ph0_smaller_50 <- r1_prob < 0.5
    # get the first index where the probability is smaller than 0.5
    first_index <- which(is_ph0_smaller_50)[1]
    print(paste("P(H0) > 0.5 for mu =", big_sim_mus[first_index]))
    lines(big_sim_mus, r1_prob, col = "black", lwd = 2)
    dev.off()
    dbDisconnect(con)
}
