setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())

#####################
##### FUNCTIONS #####
#####################

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

## derivative of bf10
der_bf10 <- function(n, mean, var1){
  (var1 * (n^2 * mean^2 * var1 + (2 * mean^2 - var1)*n-1)* 
  exp(((n * var1 / (2 * (1 + n * var1))) * zeta(n, mean)^2))) /
  (2*(1 + n * var1)^(5/2))
}

# calculates zero_points candidates of the derivative
# \frac{\sigma_1^4-2\Bar{y}^2\sigma_1^2 \pm \sqrt{ \sigma_1^4 (4\Bar{y}^4 + \sigma_1^4)}}{2\Bar{y}^2\sigma_1^4}
# mean is the mean of the sample
# var1 is the variance of the sample
# returns the two zero points of the derivative
zero_points <- function(mean, var1) {
  first_extr <- (var1 - 2 * mean^2 + sqrt(4 * mean^4 + var1^2)) /
    (2 * mean^2 * var1)
  second_extr <- (var1 - 2 * mean^2 - sqrt(4 * mean^4 + var1^2)) /
    (2 * mean^2 * var1)
  c(first_extr, second_extr)
}

# calculates the maximum points for the bayes-factor
# mean is the mean of the sample
# var1 is the variance of the sample
# returns the maximum point of the bayes factor
maxima_points <- function(mean, var){
  (mean / (sqrt((1/2)*(var + sqrt(4*mean^4 + var^2))))) *
    exp((0.5*var - mean^2 + sqrt(mean^4 + 0.25*var^2) + 
    ((2*mean^4 - mean^2 * sqrt(4*mean^4 + var^2))/var)) /
    (var + sqrt(4*mean^4 + var^2)))
}

###########################
###### CONFIGURATION ######
###########################

# define several mu's for the model m1
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
var_data <- 1
var1 <- 1
vars1 <- c(0.0001, 0.01, 1, 25, 100)
n <- seq(1, 10000, 1)
cols <- c("red", "blue", "green", "yellow", "black")


####################
##### EXAMPLES #####
####################

# print BF10 examples according to paper
# BF10 with n = 10, y_mean = 1 and var1 = [1, 25, 100, 2500, 10000]
# with bf's: [28.4, 9.2, 4.7, 0.9, 0.5]

print("BF10 Paper examples:")
tapply(vars1, print(bf10(10, 1, vars1)))
## print derivative examples
print("Derivative examples:")
tapply(mus, print(der_bf10(10, mus, var1)))
## print zero_points examples
print("Zero points examples:")
print(zero_points(1, 10000000000000000000000000000000000000000))
print(zero_points(1, vars1[1]))
## print zero_points examples
print("BF01 examples:")
print(bf10(zero_points(mus[1], var1)[1], mus[1], var1))
print("Maxima points examples (should be same as bf01):")
print(maxima_points(mus[1], var1))



###################
###### PLOTS ######
###################

# plot dependent on sample size and mu the bayes factor 01
pdf("bayes-factor-replicate.pdf")
plot(n,  bf01(bf10(n, mus[1], var1)), type = "n", xlab = "n",
     ylab = "BF_{01}", main = "Bayes Factor depending on n and mu",
     ylim = c(0, 7), xlim = c(0, 50))
legend("topright", legend = paste("mu = ", mus), col = cols, lty = 1, cex = 0.8)

for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, bf01(bf10(n, mus[i], var1)), col = cols[i], lwd = 2)
  # sanity check: are zero points correct?
  zero_point <- zero_points(mus[i], var1)[1]
  abline(v=zero_point, col = cols[i], lty = 1)
  # sanity check: can maxima points still be computed?
  #maxima_point <- bf01(bf10(zero_point, mus[i], var1))
  maxima_point <- bf01(maxima_points(mus[i], var1))
  points(zero_point, maxima_point, col = cols[i], pch = 19)
  #print(zero_points(mus[i], var1))[2]
}
dev.off()

# plot with fixed mu and variable sigma
fixed_mu <- 0.1
pdf("bf-replicate-variance.pdf")
plot(n, bf01(bf10(n,fixed_mu,vars1)), type="n", xlab="n",
ylab = "BF_{01}", main = "Bayes Factor depending on n and mu",
     ylim = c(0, 70), xlim = c(0, 1000))
legend("topright", legend = paste("var = ", vars1), col = cols, lty = 1, cex = 0.8)

for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, bf01(bf10(n, fixed_mu, vars1[i])), col = cols[i], lwd = 2)
  # sanity check: are zero points correct?
  zero_point <- zero_points(fixed_mu, vars1[i])[1]
  abline(v=zero_point, col = cols[i], lty = 1)
  # sanity check: can maxima points still be computed?
  #maxima_point <- bf01(bf10(zero_point, mus[i], var1))
  maxima_point <- bf01(maxima_points(fixed_mu, vars1[i]))
  points(zero_point, maxima_point, col = cols[i], pch = 19)
  print(zero_points(fixed_mu, vars1[i]))[2]
}
dev.off()


## plot derivatives for bf10 (mus)

pdf("bayes-factor-derivative.pdf")
plot(n,  der_bf10(n, mus[1], var1), type = "n", xlab = "n",
     ylab = "BF_{10-derivative}", main = "Bayes Factor derivative depending on n and mu",
     ylim = c(-0.5, 0.5), xlim = c(-10, 100))

legend("topright", legend = paste("mu = ", mus), col = cols, lty = 1, cex = 0.8)
abline(h = 0,  lty=2)
for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, der_bf10(n, mus[i], var1), col = cols[i], lwd = 2)
  zero_point <- zero_points(mus[i], var1)[1]
  points(zero_point, 0, col = cols[i], pch = 19)
  zero_point <- zero_points(mus[i], var1)[2]
  points(zero_point, 0, col = "grey", pch = 19)
}
dev.off()


## plot derivatives for bf10 (variances)
pdf("bayes-factor-derivative-variance.pdf")
plot(n,  der_bf10(n, fixed_mu, vars1), type = "n", xlab = "n",
     ylab = "BF_{10-derivative}", main = "Bayes Factor derivative depending on n and var",
     ylim = c(-1, 1), xlim = c(0, 10))
legend("topright", legend = paste("var = ", vars1), col = cols, lty = 1, cex = 0.8)
abline(h = 0,  lty=2)
for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, der_bf10(n, fixed_mu, vars1[i]), col = cols[i], lwd = 2)
  zero_point <- zero_points(mus[i], var1)[1]
  points(zero_point, 0, col = cols[i], pch = 19)
  zero_point <- zero_points(mus[i], var1)[2]
  points(zero_point, 0, col = "grey", pch = 19)
}
dev.off()

####################
##### LOG-PLOTS ####
####################

logscale <- c(1, 10^seq(1, 4, 0.001))

# plot dependent on sample size and mu the bayes factor 01
pdf("bayes-factor-log-scale.pdf")
plot(logscale,  bf01(bf10(logscale, mus[1], var1)), type = "n", xlab = "n",
     ylab = "BF_{01}", main = "Bayes Factor depending on n and mu",
     ylim = c(0, 7), log = "x")
legend("topright", legend = paste("mu = ", mus), col = cols, lty = 1, cex = 0.8)

for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, bf01(bf10(n, mus[i], var1)), col = cols[i], lwd = 2)
  zero_point <- zero_points(mus[i], var1)[1]
  maxima_point <- bf01(maxima_points(mus[i], var1))
  points(zero_point, maxima_point, col = cols[i], pch = 19)
}
dev.off()