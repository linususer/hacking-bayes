# Set working directory
setwd("/home/linus/git/hacking-bayes")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
library(foreach)
library(doParallel)
source("scripts/plots/define_colors.R")

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
plotProbWithBFCrit <- function(mu) {
    con <- dbConnect(duckdb(), db_file, read_only = TRUE)
    # pdf(paste("figures/realistic-sym-decision-prob-mu-", mu, ".pdf", sep = ""))
    plot(0, 0,
        xlim = c(0, 10), ylim = c(0, 1), type = "n",
        main = bquote(mu == .(mu)),
        ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(BF[crit])
    )

    r1_prob <- sapply(BF_crits, function(bf_crit) {
        as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((0.5 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND mu = ?
                                AND bf_crit = ?
                                ", list(mu, bf_crit))) / repetitions
    })

    r2_prob <- sapply(BF_crits, function(bf_crit) {
        as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND mu = ?
                                AND bf_crit = ?
                                ", list(mu, bf_crit))) / repetitions
    })
    r3_prob <- sapply(BF_crits, function(bf_crit) {
        as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((2 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND mu = ?
                                AND bf_crit = ?
                                ", list(mu, bf_crit))) / repetitions
    })

    points(BF_crits, r1_prob, col = "black", lwd = 2, pch = 4)
    text(BF_crits, r1_prob, labels = round(r1_prob, 2), pos = 3, cex = 0.8, col = "black")

    points(BF_crits, r2_prob, col = "orange", lwd = 2, pch = 4)
    text(BF_crits, r2_prob, labels = round(r2_prob, 2), pos = 3, cex = 0.8, col = "orange")

    points(BF_crits, r3_prob, col = "#009980", lwd = 2, pch = 4)
    text(BF_crits, r3_prob, labels = round(r3_prob, 2), pos = 3, cex = 0.8, col = "#009980")
    dbDisconnect(con)
    # dev.off()
}
plotProbWithBFCritBounds <- function(mu) {
    con <- dbConnect(duckdb(), db_file, read_only = TRUE)
    # Plot the decision probability for H0 given mu
    # x-axis: mu, y-axis: BF_crit
    # pdf(paste("figures/realistic-sym-decision-prob-mu-", mu, "-bf-crits.pdf", sep = ""))
    r2_prob <- sapply(BF_crits, function(bf_crit) {
        as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND mu = ?
                                AND bf_crit = ?
                                ", list(mu, bf_crit))) / repetitions
    })
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
    dbDisconnect(con)
}
# for each mu in big_sim_mus, plot the decision probability for H0 parrallel
pdf("figures/realistic-sym-decision-prob-by-bf_crit.pdf")
par(mfrow = c(2, 2))
sapply(big_sim_mus, function(mu) {
    plotProbWithBFCrit(mu)
})
legend("topright",
    legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
    col = c("black", "orange", "#009980"), lwd = 2
)
dev.off()

pdf("figures/realistic-sym-decision-prob-by-bf_crit-with-boundaries.pdf")
par(mfrow = c(2, 2))
sapply(big_sim_mus, function(mu) {
    plotProbWithBFCritBounds(mu)
})
dev.off()

big_sim_mus <- seq(0, 1, 0.01)
bf_crit <- c(3, 6, 10)
con <- dbConnect(duckdb(), db_file, read_only = TRUE)

# Plot the decision probability for H0 given mu
pdf(paste("figures/realistic-sym-decision-prob-bf-crits.pdf", sep = ""))
plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = "Decision probability for symmetrical rule in realistic setting",
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
)
bf1_prob <- sapply(big_sim_mus, function(mu) {
    dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND ROUND(mu, 6) = ROUND(?, 6)
                                AND bf_crit = ?
                                ", list(mu, bf_crit[1])) / repetitions
})
bf2_prob <- sapply(big_sim_mus, function(mu) {
    dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND ROUND(mu, 6) = ROUND(?, 6)
                                AND bf_crit = ?
                                ", list(mu, bf_crit[2])) / repetitions
})
bf3_prob <- sapply(big_sim_mus, function(mu) {
    dbGetQuery(con, "SELECT COUNT(*) FROM cauchy_sym
                                WHERE ABS((1 / SQRT(2)) - r) < 1e-6
                                AND decision = 0
                                AND ROUND(mu, 6) = ROUND(?, 6)
                                AND bf_crit = ?
                                ", list(mu, bf_crit[3])) / repetitions
})
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
dbDisconnect(con)
