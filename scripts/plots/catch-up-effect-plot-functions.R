setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
source("scripts/plots/define_colors.R")
freq_optional_stopping <- function() {
    set.seed(400)
    # draw data from normal distribution and perform two-sided t-test
    # until p-value is below 0.05 and plot p-values
    # there is no effect, but the p-value becomes significant
    data <- rnorm(2, mean = 0, sd = 1)
    p_val <- c()
    while (t.test(data)$p.value > 0.05) {
        data <- c(data, rnorm(1, mean = 0, sd = 1))
        p_val <- c(p_val, t.test(data)$p.value)
        # print(length(data))
    }
    # plot p-values
    pdf("figures/frequentistic-optional-stopping.pdf")
    plot(p_val,
        type = "l", col = "black", lwd = 4, xlab = "Sample Size", ylab = "p-value",
        main = "Frequentistic Optional Stopping",
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
    )
    points(length(p_val), t.test(data)$p.value, col = "black", pch = 4, cex = 2, lwd = 4)
    # add significance threshold with text
    text(0, 0.075, "Significance Threshold p = 0.05", pos = 4, col = "black")
    abline(h = 0.05, col = "black", lty = 2, lwd = 4)
    dev.off()
}

bayes_optional_stopping <- function() {
    # draw data from normal distribution and perform bayesian t-test
    # until the Bayes factor is higher than 3 or lower than 1/3 and plot Bayes factors
    source("bayes-factor-functions.R")
    bf_crit1 <- 1 / 3
    bf_crit2 <- 3
    data <- rnorm(1, mean = 0, sd = 1)
    bf <- bf10(length(data), mean(data), prior_var = 1)
    bf_list <- c(bf)
    while (bf > bf_crit1 && bf < bf_crit2) {
        data <- c(data, rnorm(1, mean = 0, sd = 1))
        bf <- bf10(length(data), mean(data), prior_var = 1)
        bf_list <- c(bf_list, bf)
    }
    # plot Bayes factors
    pdf("figures/bayesian-optional-stopping.pdf")
    plot(bf_list,
        type = "l", col = "#59B3E6",
        main = "Bayesian Optional Stopping",
        lwd = 2, xlab = "Sample Size", ylab = "Bayes Factor", ylim = c(1 / 3, 3)
    )
    points(length(bf_list), bf, col = "#59B3E6", pch = 4, cex = 2, lwd = 2)
    # add BF crit thresholds with text
    abline(h = 3, col = "black", lty = 2, lwd = 2)
    abline(h = 1 / 3, col = "black", lty = 2, lwd = 2)
    text(0.6, 2.9, bquote(BF[crit] * " = 3"), pos = 4, col = "black")
    text(0.6, 0.4, bquote(BF[crit] * " = 1/3"), pos = 4, col = "black")
    dev.off()
}

point_and_normal_prior <- function(true_mean = 0, observed_mean = NULL, sd_1 = 1, file_name = "figures/point-prior-and-normal-prior.pdf") {
    # plot point prior and normal prior
    pdf(file_name)
    plot(0, 0,
        xlim = c(-5, 5), ylim = c(0, 0.5), type = "n",
        main = if (is.null(observed_mean)) {
            "Point Prior vs Normal Prior"
        } else {
            bquote("Normal Prior vs Point Prior with " * bar(y) * " = " * .(observed_mean) * ", " * sigma[1] * " = " * .(sd_1))
        },
        ylab = "Density", xlab = "Mean"
    )
    lines(seq(-5, 5, 0.01), dnorm(seq(-5, 5, 0.01), true_mean, sd_1), col = "#59B3E6", lwd = 2)
    points(0, 0.5, col = "#CD1076", pch = 19, cex = 2)
    abline(v = true_mean, col = "#CD1076", lwd = 2)
    if (!is.null(observed_mean)) {
        abline(v = observed_mean, col = "#009980", lwd = 2)
        legend("topright",
            legend = c("Point Prior", "Normal Prior", bquote(bar(y) == .(observed_mean))),
            fill = c("#CD1076", "#59B3E6", "#009980")
        )
    } else {
        legend("topright",
            legend = c("Point Prior", "Normal Prior"),
            fill = c("#CD1076", "#59B3E6")
        )
    }
    dev.off()
}


bf01_mus_plot <- function(ns, mus, true_var, prior_var, logscale = FALSE, file_name = "figures/bf01-compare-mus.pdf", ylim = c(0, 10)) {
    source("bayes-factor-functions.R")
    pdf(file_name)
    str01 <- "01"
    if (logscale) {
        plot(ns, bf01(bf10_general(ns, mus[1], true_var, prior_var)),
            type = "n", xlab = "n",
            ylab = bquote(BF[.(str01)]),
            ylim = ylim,
            log = "x",
            main = bquote(BF[.(str01)] *
                " depending on n and " * bar(y) * ", " * sigma^2 * " = " * .(true_var) * ", " * sigma[1]^2 * " = " * .(prior_var))
        )
    } else {
        plot(ns, bf01(bf10_general(ns, mus[1], true_var, prior_var)),
            type = "n", xlab = "n",
            ylab = bquote(BF[.(str01)]),
            ylim = ylim,
            main = bquote(BF[.(str01)] *
                " depending on n and " * bar(y) * ", " * sigma^2 * " = " * .(true_var) * ", " * sigma[1]^2 * " = " * .(prior_var))
        )
    }
    for (i in seq_along(mus)) {
        # plot the bayes factor depending on n
        lines(ns, bf01(bf10_general(ns, mus[i], true_var, prior_var)),
            col = my_colors[i], lwd = 2
        )
    }
    legend_labels <- sapply(mus, function(mu) bquote(bar(y) == .(mu)))
    legend("topright",
        legend = legend_labels,
        col = my_colors, lty = 1, cex = 0.8
    )
    dev.off()
}

point_and_cauchy_prior <- function() {
    library("BayesFactor")
    # plot point prior and cauchy prior
    pdf("figures/point-prior-and-cauchy-prior.pdf",  width = 8, height = 4)
    # double the font size
    plot(0, 0,
        xlim = c(-3, 3), ylim = c(0, 0.5), type = "n",
        main = "Point Prior vs Cauchy Prior",
        ylab = "", xlab = "",
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
    )
    
    # Plot the Cauchy prior density
    x_vals <- seq(-5, 5, 0.01)
    lines(x_vals, dcauchy(seq(-5, 5, 0.01), 0, 1 / sqrt(2)), col = "#59B3E6", lwd = 4)
    points(0, 0.5, col = "#CD1076", pch = 19, cex = 2)
    lines(c(0, 0), y = c(0, 0.5), col = "#CD1076", lwd = 4)
    cauchy_at_0 <- dcauchy(0, 0, 1 / sqrt(2))
    points(0, cauchy_at_0, col = "black", pch = 4, cex = 2)
    # text with rscale
    #text(0, 1/sqrt(2), "r-scale", col = "black", cex = 2)
    #lines(c(-cauchy_at_0, cauchy_at_0), c(1/3, 1/3), col = "#59B3E6", lwd = 4)
    dev.off()
}


bf01_var_plot <- function(ns, mu, vars, prior_var) {
    source("bayes-factor-functions.R")
    pdf("figures/bf01-compare-vars.pdf")
    str01 <- "01"
    plot(ns, bf01(bf10_general(ns, mu, vars[1], prior_var)),
        type = "n", xlab = "n",
        ylab = bquote(BF[.(str01)]),
        ylim = c(0, 10),
        main = bquote(BF[.(str01)] *
            " depending on n and " * sigma^2 * ", " * bar(y) * " = " * .(mu) * ", " * sigma[1]^2 * " = " * .(prior_var))
    )

    for (i in seq_along(vars)) {
        # plot the bayes factor depending on n
        lines(ns, bf01(bf10_general(ns, mu, vars[i], prior_var)),
            col = my_colors[i], lwd = 2
        )
    }
    legend_labels <- sapply(vars, function(var) bquote(sigma^2 == .(var)))
    legend("topright",
        legend = legend_labels,
        col = my_colors, lty = 1, cex = 0.8
    )
    dev.off()
}

bf01_prior_var_plot <- function(ns, mu, true_var, prior_vars) {
    source("bayes-factor-functions.R")
    pdf("figures/bf01-compare-prior-vars.pdf")
    str01 <- "01"
    plot(ns, bf01(bf10_general(ns, mu, true_var, prior_vars[1])),
        type = "n", xlab = "n",
        ylab = bquote(BF[.(str01)]),
        ylim = c(0, 30),
        main = bquote(BF[.(str01)] *
            " depending on n and " * sigma[1]^2 * ", " * bar(y) * " = " * .(mu) * ", " * sigma^2 * " = " * .(true_var))
    )
    for (i in seq_along(prior_vars)) {
        # plot the bayes factor depending on n
        lines(ns, bf01(bf10_general(ns, mu, true_var, prior_vars[i])),
            col = my_colors[i], lwd = 2
        )
    }
    legend_labels <- sapply(prior_vars, function(prior_var) bquote(sigma[1]^2 == .(prior_var)))
    legend("topright",
        legend = legend_labels,
        col = my_colors, lty = 1, cex = 0.8
    )
    dev.off()
}


## plot intersection points for different BF_crits
bf01_intersection_points_bf_crit <- function(
    ns,
    mu = 0.1, true_var = 1,
    prior_var = 1, bf_crits = 3) {
    source("bayes-factor-functions.R")
    pdf("figures/bf01-intersection-points-bf-crit.pdf")
    str01 <- "01"
    plot(ns, bf01(bf10_general(ns, mu, true_var, prior_var)),
        type = "n", xlab = "n",
        ylab = bquote(BF[.(str01)]),
        main = bquote(
            "Intersection points for " * BF[.(str01)] *
                " depending on n for " * BF[crit] *
                ", " * bar(y) * " = " * .(mu) *
                ", " * sigma^2 * " = " * .(true_var) *
                ", " * sigma[1]^2 * " = " * .(prior_var)
        ), ylim = c(0, 10), xlim = c(0, 500)
    )
    lines(ns, bf01(bf10_general(ns, mu, true_var, prior_var)), col = "#009980", lwd = 2)
    for (i in seq_along(bf_crits)) {
        ns <- get_n_by_bf(mu, true_var, prior_var, bf_crits[i])
        print(ns)
        n1 <- ns[1]
        n2 <- ns[2]
        # mark n1 and n2 with a point
        points(n1, bf_crits[i], col = my_colors[i], pch = 19)
        points(n2, bf_crits[i], col = my_colors[i], pch = 19)
        # mark the area between n1 and n2 with a line
        lines(c(n1, n2), c(bf_crits[i], bf_crits[i]), col = my_colors[i], lwd = 2)
        # mark the n1 and n2 with a vertical, dashed line
        lines(c(n1, n1), c(0, bf_crits[i]), col = my_colors[i], lwd = 2, lty = 2)
        lines(c(n2, n2), c(0, bf_crits[i]), col = my_colors[i], lwd = 2, lty = 2)
    }
    legend_labels <- sapply(bf_crits, function(bf_crit) bquote(BF[crit] == .(bf_crit)))
    legend("topright",
        legend = legend_labels,
        col = my_colors, lty = 1, cex = 0.8
    )
    dev.off()
}

## plot intersection points for different mus
bf01_intersection_points_mu <- function(ns, mus, true_var = 1, prior_var = 1, bf_crit = 3) {
    source("bayes-factor-functions.R")
    pdf("figures/bf01-intersection-points-mus.pdf")
    str01 <- "01"
    plot(ns, bf01(bf10_general(ns, mus[1], true_var, prior_var)),
        type = "n", xlab = "n",
        ylab = bquote(BF[.(str01)]),
        main = bquote(
            "Intersection points for " * BF[.(str01)] *
                " depending on n for " * bar(y) *
                ", " * BF[crit] * " = " * .(bf_crit) *
                ", " * sigma^2 * " = " * .(true_var) *
                ", " * sigma[1]^2 * " = " * .(prior_var)
        ), ylim = c(0, 10), xlim = c(0, 500)
    )
    for (i in seq_along(mus)) {
        lines(ns, bf01(bf10_general(ns, mus[i], true_var, prior_var)), col = my_colors[i], lwd = 2)
        n <- get_n_by_bf(mus[i], true_var, prior_var, bf_crit)
        print(n)
        n1 <- n[1]
        n2 <- n[2]
        # mark n1 and n2 with a point
        points(n1, bf_crit, col = my_colors[i], pch = 19)
        points(n2, bf_crit, col = my_colors[i], pch = 19)
        # mark the area between n1 and n2 with a line
        lines(c(n1, n2), c(bf_crit, bf_crit), col = my_colors[i], lwd = 2)
        # mark the n1 and n2 with a vertical, dashed line
        lines(c(n1, n1), c(0, bf_crit), col = my_colors[i], lwd = 2, lty = 2)
        lines(c(n2, n2), c(0, bf_crit), col = my_colors[i], lwd = 2, lty = 2)
    }
    legend_labels <- sapply(mus, function(mu) bquote(bar(y) == .(mu)))
    legend("topright",
        legend = legend_labels,
        col = my_colors, lty = 1, cex = 0.8
    )
    dev.off()
}
