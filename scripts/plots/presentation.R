setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# import bf
source("bayes-factor-functions.R")
ns <- 1:150
fixed_mu <- 0.1
fixed_sigma <- 1

# Define the range for the normal distribution
y_norm <- seq(-5, 5, 0.01)
x_norm <- dnorm(y_norm, mean = fixed_mu, sd = 1 / sqrt(10))

# Start the plot with the functions
pdf("presentation/functions-and-normal-distribution.pdf")
plot(ns, get_y_by_bf(ns, 1, 1 / 3), type = "l", col = "skyblue", lwd = 5,
     main = "Decision Boundaries for Asymmetrical Stopping Rule",
     xlab = "Stop Count", ylab = "Estimated Mean y", xlim = c(-1, 20), ylim = c(-1, 1))
lines(ns, -get_y_by_bf(ns, 1, 1 / 3), col = "skyblue", lwd = 5)

# Find the corresponding y-values where x = 10
y_f_10 <- get_y_by_bf(10, 1, 1 / 3)
y_g_10 <- - get_y_by_bf(10, 1, 1 / 3)

abline(v = 0)
# Plot the horizontal normal distribution using x as the density and y as the range
lines(- x_norm, y_norm, col = "green", lwd = 5)

# Draw dotted lines connecting the borders
#segments(10, y_f_10, - dnorm(y_f_10, fixed_mu, 1 / sqrt(1)), y_f_10, col = "black", lty = 2)
#segments(10, y_g_10, - dnorm(y_g_10, fixed_mu, 1 / sqrt(1)), y_g_10, col = "black", lty = 2)

# Fill the area of the normal distribution
#polygon(c(0, - dnorm(y_f_10, fixed_mu, 1 / sqrt(10)),- dnorm(y_g_10, fixed_mu, 1 / sqrt(10)), 0),
 #       c(y_f_10, y_f_10, y_g_10, y_g_10), col = "green", border = NA)

abline(h = fixed_mu, col = "green", lwd = 2)
lines(ns, fixed_mu + 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 5)
lines(ns, fixed_mu - 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 5)

# Add labels and a legend
legend("topright", legend = c("Decision for H0", "Normal Distribution around mu = 0.1","95% Credibility Interval"),
       fill = c("skyblue", "green", "grey"))

dev.off()


# plot critical values for asymetric simulation
pdf("presentation/asym-critical-means.pdf")
plot(0, 0, xlim = c(-1, 20), ylim = c(-1, 1),
     main = "Decision Boundaries for Asymmetrical Stopping Rule",
     ylab = "Estimated Mean y", xlab = "Stop Count",
     type = "n")
fixed_mu <- 0.1
y_crit_h0 <- get_y_by_bf(ns, fixed_sigma, 1 / 3)
abline(v = 0)
#plot mu=0.1 with sem
lines(- x_norm, y_norm, col = "green", lwd = 5)

abline(h = fixed_mu, col = "green", lwd = 5)
lines(ns, fixed_mu + 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 5)
lines(ns, fixed_mu - 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 5)
segments(10, y_f_10, - dnorm(y_f_10, fixed_mu, 1 / sqrt(10)), y_f_10, col = "black", lty = 2)
segments(10, y_g_10, - dnorm(y_g_10, fixed_mu, 1 / sqrt(10)), y_g_10, col = "black", lty = 2)
polygon(c(0, - dnorm(y_f_10, fixed_mu, 1 / sqrt(10)),- dnorm(y_g_10, fixed_mu, 1 / sqrt(10)), 0),
        c(y_f_10, y_f_10, y_g_10, y_g_10), col = "skyblue", border = NA)
# plot h0 decision boundaries
lines(ns, y_crit_h0, col = "skyblue", lwd = 5)
lines(ns, - y_crit_h0, col = "skyblue", lwd = 5)
legend("topright", legend = c("Decision for H0", "Normal Distribution around mu = 0.1", "95% Credibility Interval"),
        fill = c("skyblue", "green", "grey"))
dev.off()

# plot critical values for asymetric simulation
pdf("presentation/asym-critical-means-1.pdf")
plot(0, 0, xlim = c(-1, 20), ylim = c(-1, 1),
     main = "Decision Boundaries for Asymmetrical Stopping Rule",
     ylab = "Estimated Mean y", xlab = "Stop Count",
     type = "n")
fixed_mu <- 0.1
y_crit_h0 <- get_y_by_bf(ns, fixed_sigma, 1 / 3)
abline(v = 0, lw = 5)
#plot mu=0.1 with sem
#abline(h = fixed_mu, col = "black", lwd = 2)
#lines(ns, fixed_mu + 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 2)
#lines(ns, fixed_mu - 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 2)
# plot h0 decision boundaries
lines(ns, y_crit_h0, col = "skyblue", lwd = 5)
lines(ns, - y_crit_h0, col = "skyblue", lwd = 5)
legend("topright", legend = c("Decision for H0"), #, "mu=0.1", "95% Credibility Interval"),
        fill = c("skyblue"))#, "black", "grey"))
dev.off()

# plot critical values for asymetric simulation
pdf("presentation/sym-critical-means.pdf")
# Find the corresponding y-values where x = 10
y_f_10 <- get_y_by_bf(10, 1, 1 / 3)
y_g_10 <- - get_y_by_bf(10, 1, 1 / 3)

y_h_10 <- get_y_by_bf(10, 1, 3)
y_i_10 <- - get_y_by_bf(10, 1, 3)

plot(0, 0, xlim = c(-1, 20), ylim = c(-1, 1),
     main = "Decision Boundaries for Symmetrical Stopping Rule",
     ylab = "Estimated Mean y", xlab = "Stop Count",
     type = "n")
fixed_mu <- 0.1
y_crit_h0 <- get_y_by_bf(ns, fixed_sigma, 1 / 3)
abline(v = 0)
#plot mu=0.1 with sem
lines(- x_norm, y_norm, col = "green", lwd = 5)
abline(h = fixed_mu, col = "green", lwd = 5)
lines(ns, fixed_mu + 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 5)
lines(ns, fixed_mu - 1.96 * fixed_sigma / sqrt(ns), col = "grey", lty = 2, lwd = 5)
# plot h0 decision boundaries
lines(ns, y_crit_h0, col = "skyblue", lwd = 5)
lines(ns, - y_crit_h0, col = "skyblue", lwd = 5)
# plot h1 decision boundaries
lines(ns, get_y_by_bf(ns, fixed_sigma, 3), col = "salmon", lwd = 5)
lines(ns, - get_y_by_bf(ns, fixed_sigma, 3), col = "salmon", lwd = 5)
# Draw dotted lines connecting the borders
segments(10, y_f_10, - dnorm(y_f_10, fixed_mu, 1 / sqrt(10)), y_f_10, col = "black", lty = 2)
segments(10, y_g_10, - dnorm(y_g_10, fixed_mu, 1 / sqrt(10)), y_g_10, col = "black", lty = 2)

# Fill the area of the normal distribution
polygon(c(0, - dnorm(y_f_10, fixed_mu, 1 / sqrt(10)),- dnorm(y_g_10, fixed_mu, 1 / sqrt(10)), 0),
        c(y_f_10, y_f_10, y_g_10, y_g_10), col = "skyblue", border = NA)

abline(v=1, col = "salmon", lwd = 5)
abline(v = 7, col = "skyblue", lwd=5)
# Draw dotted lines connecting the borders
segments(10, y_h_10, - dnorm(y_h_10, fixed_mu, 1 / sqrt(10)), y_h_10, col = "black", lty = 2)
segments(10, y_i_10, - dnorm(y_i_10, fixed_mu, 1 / sqrt(10)), y_i_10, col = "black", lty = 2)
abline(h = fixed_mu, col = "green", lwd = 5)
legend("topright", legend = c("Decision for H0", "Decision for H1", "Normal Distribution around mu = 0.1", "95% Credibility Interval"),
        fill = c("skyblue","salmon", "green", "grey"))
dev.off()