# Plots of the decision probability for H0 by BF_crit and r (thesis, symmetrical rule).
# Only reads thesis/summaries/ (written by thesis/summarise/decision-prob-bf-crit.R).
# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
source("shared/define_colors.R")
source("shared/summary-io.R")

# Initialize variables
big_sim_mus <- c(0.1, 0.2, 0.5, 0.7)

# `summary` is either the data.frame returned by probWithBFCritSummary()
# or a path to the summary .txt it wrote. Draws the panel for one mu.
plotProbWithBFCritSummary <- function(summary, mu) {
    summary_df <- read_summary(summary)
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
    summary_df <- read_summary(summary)
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

# `summary` is either the data.frame returned by decisionProbBFCritsSummary()
# or a path to the summary .txt it wrote.
plotDecisionProbBFCritsSummary <- function(summary) {
    summary_df <- read_summary(summary)
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

# for each mu in big_sim_mus, plot the decision probability for H0
bf_crit_summary <- read_summary("thesis/summaries/realistic-sym-decision-prob-by-bf_crit.txt")
pdf("thesis/figures/realistic-sym-decision-prob-by-bf_crit.pdf")
par(mfrow = c(2, 2))
invisible(sapply(big_sim_mus, function(mu) {
    plotProbWithBFCritSummary(bf_crit_summary, mu)
}))
legend("topright",
    legend = c(bquote(r == 0.5 / sqrt(2)), bquote(r == 1 / sqrt(2)), bquote(r == 2 / sqrt(2))),
    col = c("black", "orange", "#009980"), lwd = 2
)
dev.off()

pdf("thesis/figures/realistic-sym-decision-prob-by-bf_crit-with-boundaries.pdf")
par(mfrow = c(2, 2))
invisible(sapply(big_sim_mus, function(mu) {
    plotProbWithBFCritBoundsSummary(bf_crit_summary, mu)
}))
dev.off()

plotDecisionProbBFCritsSummary("thesis/summaries/realistic-sym-decision-prob-bf-crits.txt")
