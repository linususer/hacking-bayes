setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# reset seed
set.seed(NULL)
# import bf
source("bayes-factor-replicate.R")
# Rouder Simulation

mu <- 0.001

h0_count <- 0
h1_count <- 0
# run the simulation 20000 times
for (i in 1:20000) {
  h0 <- rnorm(10, mean = 0, sd = 1)
  h1 <- rnorm(10, mean = mu, sd = 1)
  bf <- bf01(bf10(10, mean(h1), var(h1)))
  # run the simulation until the bayes factor is between 1/3 and 3
  stop_count <- 1
  while (bf < 3 && bf > 1 / 3) {
    h0 <- rnorm(10, mean = 0, sd = 1)
    h1 <- rnorm(10, mean = mu, sd = 1)
    bf <- bf01(bf10(10, mean(h1), var(h1)))
    stop_count <- stop_count + 1
  }
  # count the number of times the null hypothesis is true
  if (bf > 3) {
    h0_count <- h0_count + 1
  } else {
    h1_count <- h1_count + 1
  }
  # print the results
  print(paste("Decision: ", bf, "Stop count: ", stop_count))
}
# print the results
print(paste("H0: ", h0_count))
print(paste("H1: ", h1_count))