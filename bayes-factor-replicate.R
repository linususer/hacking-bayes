setwd("/home/linus/hacking-bayes")
# clear workspace
rm(list = ls())

# calculate the bayes factor
# BF_{10} := p_k(D|M_1) / p(D|M_0)
# => BF_{10} = \frac{1}{\sqrt{1+n\sigma_1^2}} *
# e^{\frac{n\sigma_1^2}{2(1+n\sigma_1^2)}\cdot \zeta^2}

# Zeta:= \frac{\sqrt{n}(\tilde{Y}-0)}{\sigma} = \sqrt{n}\tilde{Y}
# is classic one-sample zeta test
# n is the sample size
# mean is the mean of the sample
zeta <- function(n, mean) {
  sqrt(n) * mean
}
# var1 is the variance of the sample
# n is the sample size
# mean is the mean of the sample
# bf10 is the bayes factor for the alternative hypothesis
bf10 <- function(n, mean, var1) {
  (1 / sqrt(1 + n * var1)) *
    exp(((n * var1 / (2 * (1 + n * var1))) * zeta(n, mean)^2))
}
# calculate bf01
# BF_{01} = 1 / BF_{10}
# bf10 is the bayes factor for the alternative hypothesis
# bf01 is the bayes factor for the null hypothesis
bf01 <- function(bf10) {
  1 / bf10
}
# print BF10 examples according to paper
# BF10 with n = 10, y_mean = 1 and var1 = [1, 25, 100, 2500, 10000]
# with bf's: [28.4, 9.2, 4.7, 0.9, 0.5]
vars1 <- c(1, 25, 100, 2500, 10000)
tapply(vars1, print(bf10(10, 1, vars1)))

# define several mu's for the model m1
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
var_data <- 1
var1 <- 1
n <- seq(1, 10000, 1)
cols <- c("red", "blue", "green", "yellow", "black")

# plot dependent on sample size and mu the bayes factor 01
pdf("bayes-factor-replicate.pdf")
plot(n,  bf01(bf10(n, mus[1], var1)), type = "n", xlab = "n",
     ylab = "BF_{01}", main = "Bayes Factor depending on n and mu",
     ylim = c(0, 7), xlim = c(0, 1000))
legend("topright", legend = paste("mu = ", mus), col = cols, lty = 1, cex = 0.8)

for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, bf01(bf10(n, mus[i], var1)), col = cols[i], lwd = 2)
}
dev.off()