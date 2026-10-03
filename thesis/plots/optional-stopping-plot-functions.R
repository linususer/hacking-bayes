setwd(".")
# clear workspace
rm(list = ls())
library(duckdb)
#################
### LOAD DATA ###
#################
con <- dbConnect(duckdb(), "data/hacking-bayes.duckdb")
sym <- dbGetQuery(con, "SELECT *
                        FROM bf_decision_threshold
                        WHERE symmetrical = TRUE
                        AND trial_start = 1")
asym <- dbGetQuery(con, "SELECT *
                        FROM bf_decision_threshold
                        WHERE symmetrical = FALSE
                        AND trial_start = 1")
dbDisconnect(con)
##############
### CONFIG ###
##############
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
big_sim_mus <- seq(0, 1, 0.01)
fixed_sigma <- 1
sigmas <- c(0.01, 0.1, 1, 5, 10)
max_stop_count <- 500
repetitions <- 20000

##########################
### PLOTTING FUNCTIONS ###
##########################
sim_histograms <- function(df, mus, xlim_max = 60, ylimh0 = 2500, ylimh1 = 200) {
  # plot histograms for h0 and h1 decisions
  h0_data <- df[df[, 1] == 0, , ]
  h1_data <- df[df[, 1] == 1, , ]
  df_name <- deparse(substitute(df))
  custom_breaks <- seq(0, 4000, by = 1)
  # plot both decisions
  for (mu in mus) {
    pdf(paste("thesis/figures/", df_name, "-both-", mu, ".pdf", sep = ""))
    par(mfrow = c(2, 1))
    par(mar = c(0, 5, 3, 3))
    hist(h0_data[h0_data[, 3] == mu, 2],
      main = bquote("Simulation with 20000 repetitions for " * mu == .(mu)),
      xlim = c(0, xlim_max), ylim = c(0, ylimh0),
      xlab = "",
      ylab = bquote("Decision Count (" * H[0] * ")"),
      xaxt = "n", las = 1,
      col = "#59B3E6", breaks = custom_breaks
    )
    legend("topright", legend = c(bquote(H[0]), bquote(H[1])), fill = c("#59B3E6", "#CD1076"))
    par(mar = c(5, 5, 0, 3))
    hist(h1_data[h1_data[, 3] == mu, 2],
      main = "",
      xlim = c(0, xlim_max), ylim = c(ylimh1, 0),
      xlab = "Stop Count",
      ylab = bquote("Decision Count (" * H[1] * ")"),
      las = 1,
      col = "#CD1076", breaks = custom_breaks
    )
    dev.off()
  }
}

#################################
### PLOT DECISION PROBABILITY ###
#################################
# plot the decision probability for H0 given mu and n_start
# case: "asymmetrical", "symmetrical", "both"
# n_start: number of starts for the random walk
decision_prob_curve <- function(case = "both", n_start) {
  # plot the decision probability for H0 given mu
  pdf(paste("thesis/figures/", case, if (n_start > 1) {
    paste("-n_start-", n_start, sep = "")
  }, "-decision-probability.pdf", sep = ""))
  main_text <- if (case == "both") {
    " rules"
  } else {
    " rule"
  }
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability for " * .(case) * .(main_text)),
    ylab = bquote("Decision Probability for " * H[0]), xlab = bquote(mu)
  )
  # simulation asymmetrical
  if (case == "asymmetrical" || case == "both") {
    h0_prob <- c()
    h0_data <- asym[asym[, 1] == 0, , ]
    is_ph0_smaller_50 <- FALSE
    for (mu in big_sim_mus) {
      h0_count <- length(h0_data[h0_data[, 3] == mu, 2])
      h0_prob <- c(h0_prob, h0_count / repetitions)
      if (!is_ph0_smaller_50 && h0_count / repetitions < 0.5) {
        is_ph0_smaller_50 <- TRUE
        print(paste("P(H0) > 0.5 for mu =", mu))
      }
    }
    lines(big_sim_mus, h0_prob, col = "orange", lwd = 2)
  }
  # simulation symmetrical
  if (case == "symmetrical" || case == "both") {
    h0_prob <- c()
    h0_data <- sym[sym[, 1] == 0, , ]
    is_ph0_smaller_50 <- FALSE
    for (mu in big_sim_mus) {
      h0_count <- length(h0_data[h0_data[, 3] == mu, 2])
      h0_prob <- c(h0_prob, h0_count / repetitions)
      if (!is_ph0_smaller_50 && h0_count / repetitions < 0.5) {
        is_ph0_smaller_50 <- TRUE
        print(paste("P(H0) > 0.5 for mu =", mu))
      }
    }
    lines(big_sim_mus, h0_prob, col = "#FF99FF", lwd = 2)
  }

  # legend
  if (case == "both") {
    legend("topright",
      legend = c("Asymmetrical Simulation", "Symmetrical Simulation"),
      fill = c("orange", "#FF99FF")
    )
  } else if (case == "asymmetrical") {
    legend("topright",
      legend = c("Asymmetrical Simulation"),
      fill = c("orange")
    )
  } else {
    legend("topright",
      legend = c("Symmetrical Simulation"),
      fill = c("#FF99FF")
    )
  }
  dev.off()
}
