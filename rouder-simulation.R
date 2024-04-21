setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# import bf
source("bayes-factor-functions.R")
# Rouder Simulation
fixed_mu <- 0.1
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
fixed_sigma <- 1
sigmas <- c(0.01, 0.1, 1, 5, 10)
sim <- data.frame()
h0_count <- 0
h1_count <- 0
max_stop_count <- 0

for (mu in mus){
  # reset seed
  set.seed(NULL)
  # run the simulation 20000 times
  for (i in 1:20000) {
    h0 <- rnorm(10, mean = 0, sd = 1)
    h1 <- rnorm(3, mean = mu, sd = fixed_sigma)
    bf <- bf01(bf10(length(h1), mean(h1), var(h1)))
    # run the simulation until the bayes factor is between 1/3 and 3
    stop_count <- 1
    while (bf < 3 && bf > 1 / 3) {
      h0 <- rnorm(10, mean = 0, sd = 1)
      h1 <- c(h1, rnorm(1, mean = mu, sd = 1))
      bf <- bf01(bf10(length(h1), mean(h1), var(h1)))
      stop_count <- stop_count + 1
    }
    # count the number of times the null hypothesis is true
    # and save the results for plotting
    if (bf > 3) {
      h0_count <- h0_count + 1
      sim <- rbind(sim, c(0, stop_count, mu))
    } else {
      h1_count <- h1_count + 1
      sim <- rbind(sim, c(1, stop_count, mu))
    }
    # determine new maximum stop count
    if (stop_count > max_stop_count) {
      max_stop_count <- stop_count
    }
    # print the results
    # print(paste("Decision: ", bf, "Stop count: ", stop_count))
  }
  # print the results
  print(paste("Simulation is done for ", mu))
  print(paste("H0: ", h0_count))
  print(paste("H1: ", h1_count))
}
################
### PLOTTING ###
################

h0_data <- sim[sim[, 1] == 0, 2:3]
h1_data <- sim[sim[, 1] == 1, 2:3]
hist_data_h0 <- hist(h0_data[h0_data[, 2] == fixed_mu, 1], breaks = 50,
                     plot = FALSE)
hist_data_h1 <- hist(h1_data[h1_data[, 2] == fixed_mu, 1], breaks = 50,
                     plot = FALSE)
# plot only h0 decisions
pdf("plots/rs-h0.pdf")
hist(h0_data[h0_data[, 2] == fixed_mu, 1], col = "skyblue", breaks = 50,
     xlim = c(0, 150), ylim = c(0, 5000), main = "Rouder Simulation for h0",
     xlab = "Stop Count", ylab = "Decision Count")
dev.off()
#plot only h1 decisions
pdf("plots/rs-h1.pdf")
hist(h1_data[h1_data[, 2] == fixed_mu, 1], col = "salmon", breaks = 50,
     xlim = c(0, 150), ylim = c(0, 5000), main = "Rouder Simulation for h1",
     xlab = "Stop Count", ylab = "Decision Count")
dev.off()
# plot both decisions
pdf("plots/rs-both.pdf")
for (mu in mus) {
  hist_data_h1_mirr <- hist(h1_data[h1_data[, 2] == mu, 1],
                            breaks = 50, plot = FALSE)
  hist_data_h1_mirr$counts <- -hist_data_h1$counts
  plot(hist_data_h1_mirr, xlim = c(0, 150), ylim = c(-5000, 5000),
       col = "salmon",
       main = "Rouder Simulation for h0 vs h1",
       xlab = "Stop Count",
       ylab = "Decision Count (h0 vs. h1)")
  hist(h0_data[h0_data[, 2] == mu, 1], col = "skyblue",
       breaks = 50, add = TRUE)
  legend("topright", legend = c("h0", "h1"), fill = c("skyblue", "salmon"))
}
dev.off()
# plot both decisions as a relative frequency
pdf("plots/rs-both-rel.pdf")
for (mu in mus) {
  plot(0, 0, xlim = c(0, 150), ylim = c(0, 100), type = "n",
       xlab = "Stop Count", ylab = "Relative Decision Count",
       main = paste("Rouder Simulation for h0 vs h1; mu=", mu))
  h_rel_arr <- c()
  for (i in 1:max_stop_count) {
    h0_count <- length(h0_data[h0_data[, 1] == i & h0_data[, 2] == mu, 1])
    h1_count <- length(h1_data[h1_data[, 1] == i & h1_data[, 2] == mu, 1])
    h_rel <- h0_count / h1_count
    h_rel_arr <- c(h_rel_arr, h_rel)
    points(i, h0_count, col = "skyblue")
    points(i, h1_count, col = "salmon")
    points(i, h_rel, col = "black")
  }
  lines(h_rel_arr, col = "black")
  legend("topright", legend = c("h0", "h1", "h0/h1"),
         fill = c("skyblue", "salmon", "black"))
}
dev.off()

# plot mu regarding the probability p(H_0),
# so the decision probability for H_0.

pdf("plots/rs-decision-probability-h0.pdf")
plot(0, 0, xlim = c(0,1), ylim = c(0,1), type = "n")
for (mu in mus) {
  h0_count <- length(h0_data[h0_data[, 1] == i && h0_data[, 2] == mu, 1])
  h1_count <- length(h1_data[h1_data[, 1] == i && h1_data[, 2] == mu, 1])
  h0_prob <- h0_count / (h1_count + h0_count)
  lines(mu, h0_prob)
}

# plot critical mean values to estimate decision boundaries 
# dependent on effect size
pdf("plots/rs-critical-means.pdf")
plot(1:1000, get_y_by_bf(1:1000, fixed_sigma, 3)[2])
# lines(0:1000, get_y_by_bf(0:1000, fixed_sigma, 3)[2])
dev.off()