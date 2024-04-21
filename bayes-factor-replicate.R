setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# import functions
source("bayes-factor-functions.R")

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
bf <- 1/3


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
print("Newtons Method")
# print(get_n_by_bf(0.1, 1, 3))


###################
###### PLOTS ######
###################

# plot dependent on sample size and mu the bayes factor 01
pdf("plots/bf.pdf")
plot(n,  bf01(bf10(n, mus[1], var1)), type = "n", xlab = "n",
     ylab = "BF_{01}", main = "Bayes Factor depending on n and mu",
     ylim = c(0, 7), xlim = c(0, 200))
legend("topright", legend = paste("mu = ", mus), col = cols, lty = 1, cex = 0.8)

for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, bf01(bf10(n, mus[i], var1)), col = cols[i], lwd = 2)
  # sanity check: are zero points correct?
  zero_point <- zero_points(mus[i], var1)[1]
  #abline(v=zero_point, col = cols[i], lty = 1)
  # sanity check: can maxima points still be computed?
  #maxima_point <- bf01(bf10(zero_point, mus[i], var1))
  maxima_point <- bf01(maxima_points(mus[i], var1))
  points(zero_point, maxima_point, col = cols[i], pch = 19)
  #print(zero_points(mus[i], var1))[2]
  n_tuple <- get_n_by_bf(mus[i], var1, bf)
  points(n_tuple, rep(bf, 2), pch = 19, col = cols[i])
}
abline(h = bf)
dev.off()

# plot with fixed mu and variable sigma
fixed_mu <- 0.1
pdf("plots/bf-variance.pdf")
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

pdf("plots/bf-derivative.pdf")
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
pdf("plots/bf-derivative-variance.pdf")
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
pdf("plots/bf-log-scale.pdf")
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

# log ( bf ) 
print("Log Bayes Factor examples:")
tapply(mus, print(log(bf10(10, mus, var1))))
pdf("bf-log.pdf")
plot(n, log(bf10(n, mus[1], var1)), type = "n", xlab = "n",
     ylab = "log(BF_{10})", main = "log(Bayes Factor depending on n and mu",
     xlim = c(0, 50),
     ylim = c(-2, 2))
legend("topright", legend = paste("mu = ", mus), col = cols, lty = 1, cex = 0.8)
for (i in 1:5) {
  # plot the bayes factor depending on n
  lines(n, log(bf01(bf10(n, mus[i], var1))), col = cols[i], lwd = 2)
  zero_point <- zero_points(mus[i], var1)[1]
  maxima_point <- bf01(maxima_points(mus[i], var1))
  abline(v=zero_point, col = cols[i], lty = 1)
  points(zero_point, maxima_point, col = cols[i], pch = 19)
}
dev.off()