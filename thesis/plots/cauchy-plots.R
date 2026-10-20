# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
library(foreach)
library(doParallel)
source("shared/define_colors.R")

# Register parallel backend
registerDoParallel(cores = detectCores())

# Initialize variables
big_sim_mus <- c(0.1, 0.2, 0.5, 0.7)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

# Convert CSV to SQLite DB
db_file <- "data/hacking-bayes.duckdb"

# Load results # does not work because file is too big
# sym <- read.csv("data/realistic_sym_simulation.csv")
# Plot the decision probability for H0 given mu
# x-axis: BF_crit, y-axis: P(H0), color: r
# Summary: decision probability for H0 per mu, r and bf_crit.
# Used by plotProbWithBFCritSummary() and plotProbWithBFCritBoundsSummary().
probWithBFCritSummary <- function(mus, summary_path = NULL) {
    con <- dbConnect(duckdb(), db_file, read_only = TRUE)
    summary_rows <- list()
    for (mu in mus) {
        r1_prob <- sapply(BF_crits, function(bf_crit) {
            as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((0.5 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND mu = ?
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit))) / repetitions
        })

        r2_prob <- sapply(BF_crits, function(bf_crit) {
            as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND mu = ?
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit))) / repetitions
        })
        r3_prob <- sapply(BF_crits, function(bf_crit) {
            as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((2 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND mu = ?
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit))) / repetitions
        })
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            mu = mu,
            r = rep(c(0.5, 1, 2) / sqrt(2), each = length(BF_crits)),
            bf_crit = rep(BF_crits, 3),
            prob = c(r1_prob, r2_prob, r3_prob)
        )
    }
    dbDisconnect(con)
    summary_df <- do.call(rbind, summary_rows)
    rownames(summary_df) <- NULL

    if (!is.null(summary_path)) {
        dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
        write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
        cat("=> Saved bf_crit summary to:", summary_path, "\n")
    }

    summary_df
}

# `summary` is either the data.frame returned by probWithBFCritSummary()
# or a path to the summary .txt it wrote. Draws the panel for one mu.
plotProbWithBFCritSummary <- function(summary, mu) {
    summary_df <- if (is.character(summary)) {
        read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
        summary
    }
    get_prob <- function(r) {
        summary_df[abs(summary_df$mu - mu) < 1e-9 & abs(summary_df$r - r) < 1e-6, ]
    }
    # pdf(paste("thesis/figures/realistic-sym-decision-prob-mu-", mu, ".pdf", sep = ""))
    plot(0, 0,
        xlim = c(0, 10), ylim = c(0, 1), type = "n",
        main = bquote(mu == .(mu)),
        ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(BF[crit])
    )

    r1 <- get_prob(0.5 / sqrt(2))
    r2 <- get_prob(1 / sqrt(2))
    r3 <- get_prob(2 / sqrt(2))
    BF_crits <- r1$bf_crit
    r1_prob <- r1$prob
    r2_prob <- r2$prob
    r3_prob <- r3$prob

    points(BF_crits, r1_prob, col = "black", lwd = 2, pch = 4)
    text(BF_crits, r1_prob, labels = round(r1_prob, 2), pos = 3, cex = 0.8, col = "black")

    points(BF_crits, r2_prob, col = "orange", lwd = 2, pch = 4)
    text(BF_crits, r2_prob, labels = round(r2_prob, 2), pos = 3, cex = 0.8, col = "orange")

    points(BF_crits, r3_prob, col = "#009980", lwd = 2, pch = 4)
    text(BF_crits, r3_prob, labels = round(r3_prob, 2), pos = 3, cex = 0.8, col = "#009980")
    # dev.off()
    invisible(summary_df)
}
# `summary` is either the data.frame returned by probWithBFCritSummary()
# or a path to the summary .txt it wrote. Only uses r = 1 / sqrt(2).
plotProbWithBFCritBoundsSummary <- function(summary, mu) {
    summary_df <- if (is.character(summary)) {
        read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
        summary
    }
    # Plot the decision probability for H0 given mu
    # x-axis: mu, y-axis: BF_crit
    # pdf(paste("thesis/figures/realistic-sym-decision-prob-mu-", mu, "-bf-crits.pdf", sep = ""))
    r2 <- summary_df[abs(summary_df$mu - mu) < 1e-9 & abs(summary_df$r - 1 / sqrt(2)) < 1e-6, ]
    BF_crits <- r2$bf_crit
    r2_prob <- r2$prob
    plot(0, 0,
        xlim = c(0, 10), ylim = c(0, 1), type = "n",
        main = bquote(mu == .(mu)),
        xlab = bquote(BF[crit]), ylab = bquote(P(H[0]))
    )
    points(BF_crits, r2_prob, col = "black", lwd = 2, pch = 4)
    text(BF_crits, r2_prob, labels = round(r2_prob, 2), pos = 3, cex = 0.8, col = "black")
    lines(1:10, sapply(1:10, function(x) {
        1 / x
    }), col = "black", lty = 2)
    # dev.off()
    invisible(summary_df)
}
# Summary: decision probability for H0 per mu and bf_crit (r = 1 / sqrt(2)).
decisionProbBFCritsSummary <- function(big_sim_mus, bf_crit, summary_path = NULL) {
    con <- dbConnect(duckdb(), db_file, read_only = TRUE)
    bf1_prob <- sapply(big_sim_mus, function(mu) {
        dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND ROUND(mu, 6) = ROUND(?, 6)
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit[1])) / repetitions
    })
    bf2_prob <- sapply(big_sim_mus, function(mu) {
        dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND ROUND(mu, 6) = ROUND(?, 6)
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit[2])) / repetitions
    })
    bf3_prob <- sapply(big_sim_mus, function(mu) {
        dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                    WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                    AND decision = 0
                                    AND ROUND(mu, 6) = ROUND(?, 6)
                                    AND bf_crit = ?
                                    AND trial_start = 2 AND trial_end = 100000
                                    ", list(mu, bf_crit[3])) / repetitions
    })
    dbDisconnect(con)
    summary_df <- data.frame(
        mu = rep(big_sim_mus, 3),
        bf_crit = rep(bf_crit, each = length(big_sim_mus)),
        prob = c(as.numeric(unlist(bf1_prob)), as.numeric(unlist(bf2_prob)), as.numeric(unlist(bf3_prob)))
    )

    if (!is.null(summary_path)) {
        dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
        write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
        cat("=> Saved bf_crits summary to:", summary_path, "\n")
    }

    summary_df
}

# `summary` is either the data.frame returned by decisionProbBFCritsSummary()
# or a path to the summary .txt it wrote.
plotDecisionProbBFCritsSummary <- function(summary) {
    summary_df <- if (is.character(summary)) {
        read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
        summary
    }
    bf_crit <- unique(summary_df$bf_crit)
    big_sim_mus <- summary_df$mu[summary_df$bf_crit == bf_crit[1]]
    bf1_prob <- summary_df$prob[summary_df$bf_crit == bf_crit[1]]
    bf2_prob <- summary_df$prob[summary_df$bf_crit == bf_crit[2]]
    bf3_prob <- summary_df$prob[summary_df$bf_crit == bf_crit[3]]
    # Plot the decision probability for H0 given mu
    pdf(paste("thesis/figures/realistic-sym-decision-prob-bf-crits.pdf", sep = ""))
    plot(0, 0,
        xlim = c(0, 1), ylim = c(0, 1), type = "n",
        main = "Decision probability for symmetrical rule in realistic setting",
        ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
    )
    # lines(big_sim_mus, r1_prob, col = "black", lwd = 2)
    lines(big_sim_mus, bf1_prob, col = "black", lwd = 2)
    lines(big_sim_mus, bf2_prob, col = "orange", lwd = 2)
    lines(big_sim_mus, bf3_prob, col = "#009980", lwd = 2)
    # boundaries
    abline(h = 1 / bf_crit[1], col = "black", lty = 2)
    abline(h = 1 / bf_crit[2], col = "orange", lty = 2)
    abline(h = 1 / bf_crit[3], col = "#009980", lty = 2)

    # lines(big_sim_mus, r3_prob, col = "#009980", lwd = 2)
    legend("topright",
        legend = c(bquote(BF[crit] == 3), bquote(BF[crit] == 6), bquote(BF[crit] == 10)),
        col = c("black", "orange", "#009980"), lwd = 2
    )
    dev.off()
    invisible(summary_df)
}

# for each mu in big_sim_mus, plot the decision probability for H0 parrallel
bf_crit_summary <- probWithBFCritSummary(
    big_sim_mus,
    summary_path = "thesis/summaries/realistic-sym-decision-prob-by-bf_crit.txt"
)
pdf("thesis/figures/realistic-sym-decision-prob-by-bf_crit.pdf")
par(mfrow = c(2, 2))
sapply(big_sim_mus, function(mu) {
    plotProbWithBFCritSummary(bf_crit_summary, mu)
})
legend("topright",
    legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
    col = c("black", "orange", "#009980"), lwd = 2
)
dev.off()

pdf("thesis/figures/realistic-sym-decision-prob-by-bf_crit-with-boundaries.pdf")
par(mfrow = c(2, 2))
sapply(big_sim_mus, function(mu) {
    plotProbWithBFCritBoundsSummary(bf_crit_summary, mu)
})
dev.off()

big_sim_mus <- seq(0, 1, 0.01)
bf_crit <- c(3, 6, 10)
bf_crits_summary <- decisionProbBFCritsSummary(
    big_sim_mus, bf_crit,
    summary_path = "thesis/summaries/realistic-sym-decision-prob-bf-crits.txt"
)
plotDecisionProbBFCritsSummary(bf_crits_summary)
