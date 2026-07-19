library(safestats)

set.seed(5224)

# Standardbeispiel, das auch in 2. der Vignette verwendet wird.



## 1. Design a one-sample safe t-test (this fixes the e-value's shape,
##    i.e. its 'deltaS' parameter, and the decision threshold 1/alpha)
alpha    <- 0.05
deltaMin <- 0.1  # minimal effect size we want to be able to detect

designObj <- designSafeT(
  deltaMin    = deltaMin,
  alpha       = alpha,
  beta        = 0.2,
  alternative = "greater",
  testType    = "oneSample"
)
designObj

nMax <- 2000   # max sample size

## 2. Simulate data one observation at a time (here under H0: mean = 0)
##    and track the e-value as each new data point comes in --
##    this sequence *is* the random walk / test martingale.
x        <- rnorm(nMax, mean = 0, sd = 1)
eValues  <- rep(NA_real_, nMax)
stopped  <- FALSE
stopN    <- NA

for (n in 2:nMax) {                       # t-stat needs at least 2 obs
  tStat <- t.test(x[1:n])$statistic

  eValues[n] <- safeTTestStat(
    t           = tStat,
    parameter   = designObj$parameter,    # deltaS from the design
    n1          = n,
    alternative = "greater"
  )

  if (!stopped && eValues[n] > 1 / alpha) {
    stopped <- TRUE
    stopN   <- n
  }
}

## 3. Plot the e-value random walk with the optional-stopping boundary
plot(2:nMax, eValues[2:nMax], type = "l", lwd = 2,
     xlab = "n (sample size)", ylab = "e-value",
     main = "E-value process under optional stopping")
abline(h = 1 / alpha, col = "red", lty = 2, lwd = 2)

if (stopped) {
  points(stopN, eValues[stopN], pch = 19, col = "red")
  legend("topleft",
    legend = c(paste0("1/alpha = ", 1/alpha),
               paste0("stopped at n = ", stopN)),
    col = c("red", "red"), lty = c(2, NA), pch = c(NA, 19), bty = "n")
} else {
  legend("topleft", legend = paste0("1/alpha = ", 1/alpha),
         col = "red", lty = 2, bty = "n")
}

pdf("figures/report/e-value-random-walk.pdf", width = 8, height = 6)

plot(2:nMax, eValues[2:nMax], type = "l", lwd = 2,
     xlab = "n (sample size)", ylab = "e-value",
     main = "E-value process under optional stopping")
abline(h = 1 / alpha, col = "red", lty = 2, lwd = 2)

if (stopped) {
  points(stopN, eValues[stopN], pch = 19, col = "red")
  legend("topleft",
    legend = c(paste0("1/alpha = ", 1/alpha),
               paste0("stopped at n = ", stopN)),
    col = c("red", "red"), lty = c(2, NA), pch = c(NA, 19), bty = "n")
} else {
  legend("topleft", legend = paste0("1/alpha = ", 1/alpha),
         col = "red", lty = 2, bty = "n")
}

dev.off()