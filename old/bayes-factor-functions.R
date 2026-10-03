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
# true_var is the known variance of the true distribution
# bf10 is the bayes factor for the alternative hypothesis
# in comparison to the null hypothesis
bf10_general <- function(n, mean, true_var, prior_var) {
  (sqrt(true_var) / sqrt(true_var + n * prior_var)) *
    exp(((n * prior_var / (2 * true_var * (true_var + n * prior_var))) * zeta(n, mean)^2))
}

# var1 is the variance of the sample
# n is the sample size
# mean is the mean of the sample
# bf10 is the bayes factor for the alternative hypothesis with known variance 1
# in comparison to the null hypothesis
bf10 <- function(n, mean, prior_var) {
  bf10_general(n, mean, true_var = 1, prior_var)
}

# calculate bf01
# BF_{01} = 1 / BF_{10}
# bf10 is the bayes factor for the alternative hypothesis
# bf01 is the bayes factor for the null hypothesis
bf01 <- function(bf10) {
  1 / bf10
}

## derivative of bf10
der_bf10 <- function(n, mean, var1) {
  (var1 * (n^2 * mean^2 * var1 + (2 * mean^2 - var1) * n - 1) *
    exp(((n * var1 / (2 * (1 + n * var1))) * zeta(n, mean)^2))) /
    (2 * (1 + n * var1)^(5 / 2))
}
der_bf10_general <- function(n, mean, var, prior_var) {
  (prior_var * (n^2 * mean^2 * prior_var + (2 * mean^2 * var - var * prior_var) * n - var^2) *
    exp(((n * prior_var / (2 * var * (var + n * prior_var))) * zeta(n, mean)^2))) /
    (2 * var * (var + n * prior_var)^(5 / 2))
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

zero_points_general <- function(mean, var, prior_var) {
  first_extr <- var * (prior_var - 2 * mean^2 + sqrt(4 * mean^4 + prior_var^2)) /
    (2 * mean^2 * prior_var)
  second_extr <- var * (prior_var - 2 * mean^2 - sqrt(4 * mean^4 + prior_var^2)) /
    (2 * mean^2 * prior_var)
  c(first_extr, second_extr)
}

# calculates the maximum points for the bayes-factor
# mean is the mean of the sample
# var1 is the variance of the sample
# returns the maximum point of the bayes factor
maxima_points <- function(mean, var) {
  (mean / (sqrt((1 / 2) * (var + sqrt(4 * mean^4 + var^2))))) *
    exp((0.5 * var - mean^2 + sqrt(mean^4 + 0.25 * var^2) +
      ((2 * mean^4 - mean^2 * sqrt(4 * mean^4 + var^2)) / var)) /
      (var + sqrt(4 * mean^4 + var^2)))
}

# calculate n with given BF and newtons method
# mean is the mean of the sample
# var is the variance of the sample
# bf is the wanted Bayes-Factor
# returns vector of c(n_1,n_2) or c(NaN, NaN)
get_n_by_bf <- function(mean, var, prior_var, bf) {
  extrema <- zero_points_general(mean, var, prior_var)[1]
  n_1 <- ceiling(extrema - 0.001)
  n_2 <- floor(extrema + 0.001)

  if (bf > bf01(bf10_general(n_1, mean, var, prior_var)) ||
    bf > bf01(bf10_general(n_2, mean, var, prior_var))) {
    return(c(NaN, NaN))
  }
  # calculate left boundary
  while (bf < bf01(bf10_general(n_1, mean, var, prior_var)) && n_1 >= 1) {
    # print(n_1)
    n_1 <- n_1 + bf01(bf10_general(n_1, mean, var, prior_var)) / bf01(der_bf10_general(n_1, mean, var, prior_var))
    # if (n_1 < 1) {
    #  n_1 <- NaN
    # }
  }
  # calculate right boundary
  while (bf < bf01(bf10_general(n_2, mean, var, prior_var))) {
    # print(n_2)
    n_2 <- n_2 + bf01(bf10_general(n_2, mean, var, prior_var)) / bf01(der_bf10_general(n_2, mean, var, prior_var))
  }
  c(n_1, n_2)
}

# calculate \Bar{y}_{crit} with given BF, var and n
# n is the sample size
# prior_var is the alternative prior variance
# bf is the Bayes-Factor threshold
# return vector of +y_crit

get_y_by_bf <- function(n, prior_var, bf) {
  y <- sqrt((2 * log(bf) + log(1 + n * prior_var)) / (n^2 * prior_var)
    + (2 * log(bf) + log(1 + n * prior_var)) / (n))
  y
}


# error phi function
# n is the stop count
# mean is the true mean of the sample
# var is the variance of the mean
# ycrit is the critical value
# returns the error probability
phi <- function(n, mean, var, ycrit) {
  integrand <- function(x) {
    exp(-(x - mean)^2 / (2 * var))
  }
  # print(cat("Integral with ycrit: ", ycrit, " and n: ", n))
  integral <- integrate(f = integrand, lower = -Inf, upper = ycrit)
  (1 / sqrt(2 * pi * var)) * integral$value
}

# P(H0 | N = n)
# n is the stop count
# mean is the true mean of the sample
# var is the variance of the mean
# bf_crit is the Bayes-Factor threshold for BF_10!
# returns the stopping probability at n
h0_stopping_probability_at_n <- function(n, mean, prior_var, bf_crit) {
  var <- (1 / sqrt(n))^2
  y_crit <- get_y_by_bf(n, prior_var, bf_crit)
  if (is.nan(y_crit)) {
    0
  } else {
    phi(n, mean, var, y_crit) -
      phi(n, mean, var, -y_crit)
  }
}

# calculate indecisive probability until stop count k
# mean is the true mean of the sample
# var is the variance of the mean
# bf_crit is the Bayes-Factor threshold for BF_10!
# returns the indecisive probability until stop count k
asym_indecisive_probability_until <- function(k, mean, prior_var, bf_crit, start) {
  product <- 1
  i <- start
  while (i <= k) {
    product <- product * (1 - h0_stopping_probability_at_n(
      i, mean,
      prior_var, bf_crit
    ))
    i <- i + 1
  }
  product
}

# calculate the decision probability for H0
# k is the current stop count
# mean is the true mean of the sample
# bf_crit is the Bayes-Factor threshold for BF_10!
# return list with the probabilities for H0 until stop count k
asym_decision_probability <-
  function(mean, prior_var, bf_crit, start, verbose = FALSE) {
    bf_crit <- 1 / bf_crit # H0
    lst <- vector()
    sum <- 0
    sum_diff <- 1
    k <- start
    while (k <= 150) { # && (k <= 15 ||sum_diff > 10^-20)) {
      sum_diff <- h0_stopping_probability_at_n(k, mean, prior_var, bf_crit) *
        asym_indecisive_probability_until(k - 1, mean, prior_var, bf_crit, start)
      sum <- sum + sum_diff
      # add to data frame
      lst[k] <- sum
      if (verbose) {
        print(paste("Summe:", sum, "k:", k, "sum_diff:", sum_diff))
      }
      k <- k + 1
    }
    lst
  }

# a <- asym_decision_probability(0.1, 1, 3)
# print(a)

##########################
##### SYMMETRIC CASE #####
##########################

h1_stopping_probability_at_n <- function(n, mean, prior_var, bf_crit) {
  var <- (1 / sqrt(n))^2
  y_crit <- get_y_by_bf(n, prior_var, bf_crit)
  if (is.nan(y_crit)) {
    0
  } else {
    1 - phi(n, mean, var, y_crit) +
      phi(n, mean, var, -y_crit)
  }
}

sym_indecisive_probability_until <- function(k, mean, prior_var,
                                             bf_crit1, bf_crit2, start) {
  product <- 1
  i <- start
  # while (h0_stopping_probability_at_n(i, mean, prior_var, bf_crit1) == 0) {
  #   i <- i + 1
  # }
  while (i <= k) {
    product <- product *
      (1 - h0_stopping_probability_at_n(i, mean, prior_var, bf_crit1) -
        h1_stopping_probability_at_n(i, mean, prior_var, bf_crit2))
    i <- i + 1
  }
  product
}
# calculate the decision probability for H0
sym_decision_probability <-
  function(mean, prior_var, bf_crit, start, verbose = FALSE) {
    bf_crit1 <- 1 / bf_crit # H0
    bf_crit2 <- bf_crit # H1
    lst <- vector()
    sum <- 0
    sum_diff <- 1
    k <- start
    # while (h0_stopping_probability_at_n(k, mean, prior_var, bf_crit1) == 0) {
    #   k <- k + 1
    # }
    while (k <= 150) { # && (k <= 15 || sum_diff > 10^-20)) {
      sum_diff <- h0_stopping_probability_at_n(k, mean, prior_var, bf_crit1) *
        sym_indecisive_probability_until(
          k - 1, mean, prior_var,
          bf_crit1, bf_crit2, start
        )
      sum <- sum + sum_diff
      # add to data frame
      lst[k] <- sum
      if (verbose) {
        print(paste("Summe:", sum, "k:", k, "sum_diff:", sum_diff))
      }
      k <- k + 1
    }
    lst
  }
# b <- sym_decision_probability(0.1, 1, 3)
# print(b)
