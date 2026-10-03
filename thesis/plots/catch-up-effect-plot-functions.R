setwd(".")
# clear workspace
rm(list = ls())
source("shared/define_colors.R")
freq_optional_stopping <- function() {
    # set.seed(400)
    # draw data from normal distribution and perform two-sided t-test
    # until p-value is below 0.05 and plot p-values
    # there is no effect, but the p-value becomes significant
    data <- rnorm(2, mean = 0, sd = 1)
    p_val <- c()
    count <- 2
    while (count <= 2000) {
        data <- c(data, rnorm(1, mean = 0, sd = 1))
        p_val <- c(p_val, t.test(data)$p.value)
        # print(length(data))
        count <- count + 1
    }
    # plot p-values
    pdf("thesis/figures/frequentistic-optional-stopping.pdf")
    plot(p_val,
        type = "l", col = "black", lwd = 4, xlab = "Sample Size", ylab = "p-value",
        main = "Standard Frequentist t-test", xlim = c(0,2000), ylim = c(0, 1),
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
    )
    # points(length(p_val), t.test(data)$p.value, col = "black", pch = 4, cex = 2, lwd = 4)
    # add significance threshold with text
    text(1170, 0.075, "Significance threshold p = 0.05", pos = 4, col = "black")
    abline(h = 0.05, col = "black", lty = 2, lwd = 4)
    dev.off()
}

point_and_cauchy_prior <- function() {
    library("BayesFactor")
    # plot point prior and cauchy prior
    pdf("thesis/figures/point-prior-and-cauchy-prior.pdf", width = 8, height = 4)
    # double the font size
    plot(0, 0,
        xlim = c(-3, 3), ylim = c(0, 0.5), type = "n",
        # main = "Point Prior vs Cauchy Prior",
        ylab = "", xlab = bquote("Effect Size " * delta),
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
    )

    # Plot the Cauchy prior density
    x_vals <- seq(-5, 5, 0.01)
    lines(x_vals, dcauchy(seq(-5, 5, 0.01), 0, 1 / sqrt(2)), col = "#59B3E6", lwd = 4)
    points(0, 0.5, col = "#CD1076", pch = 19, cex = 2)
    lines(c(0, 0), y = c(0, 0.5), col = "#CD1076", lwd = 4)
    cauchy_at_0 <- dcauchy(0, 0, 1 / sqrt(2))
    # points(0, cauchy_at_0, col = "black", pch = 4, cex = 2) #
    arrows(x0 = 0, x1 = 1 / sqrt(2), y0 = dcauchy(1 / sqrt(2)), y1 = dcauchy(1 / sqrt(2)), length = 0.15, lwd = 4)
    text(x = 1 / (2 * sqrt(2)) - 0.05, y = dcauchy(1 / sqrt(2)) + 0.03, label = "r")
    legend("topleft", legend = c(bquote(H[0]:delta == 0), bquote(H[1]:delta %~% Cauchy(r))), fill = c("#CD1076", "#59B3E6"))
    text(x = 2.5, y = 0.45, labels = bquote(BF["01"] == frac(p(data * " | " * H[0]), p(data * " | " * H[1]))))
    dev.off()
}
