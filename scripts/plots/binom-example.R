# Example of two added binomial distributions with p=0.25 and p=0.75 plotted
pdf("figures/binom-example.pdf")
plot(
    x = 0:20,
    y = dbinom(0:20, 20, 0.25),
    type = "h",
    col = "blue",
    lwd = 2,
    xlab = "Number of heads",
    ylab = "Probability",
    main = "Binomial distribution with p=0.25 and p=0.75"
)
lines(
    x = 0:20,
    y = dbinom(0:20, 20, 0.75),
    type = "h",
    col = "red",
    lwd = 2
)
legend(
    "topright",
    legend = c("p=0.25", "p=0.75"),
    col = c("blue", "red"),
    lwd = 2
)
dev.off()

# Example of two added binomial distributions with p=0.25 and p=0.75 plotted but only weighted with 0.5

pdf("figures/binom-example-weighted.pdf")
plot(
    x = 0:20,
    y = 0.5 * dbinom(0:20, 20, 0.25) + 0.5 * dbinom(0:20, 20, 0.75),
    type = "h",
    col = "black",
    lwd = 2,
    xlab = "Number of heads",
    ylab = "Probability",
    main = "Weighted binomial distribution with p=0.25 and p=0.75"
)
dev.off()

# Weighted example with 0.25 and 0.75 in comparison to 0.5
pdf("figures/binom-example-weighted-comp.pdf")
plot(
    x = 0:20,
    y = dbinom(0:20, 20, 0.5),
    type = "h",
    col = "black",
    lwd = 2,
    xlab = "Number of heads",
    ylab = "Probability",
    main = "Binomial distribution with p=0.5"
)
lines(
    x = 0:20,
    y = 0.5 * dbinom(0:20, 20, 0.25) + 0.5 * dbinom(0:20, 20, 0.75),
    type = "h",
    col = "red",
    lwd = 2
)
legend(
    "topright",
    legend = c("p=0.5", "p=0.25, p=0.75"),
    col = c("black", "red"),
    lwd = 2
)
dev.off()
