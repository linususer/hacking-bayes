setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
source("bayes-factor-functions.R")
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
    pdf(paste("figures/", df_name, "-both-", mu, ".pdf", sep = ""))
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
################
### PLOTTING ###
################

# plot critical mean values to estimate decision boundaries
# dependent on effect size
# sample size ns
# prior variance prior_var
# bf_crit Bayes factor threshold for H1
# step actual step for the normal distribution
# symmetric if the Optional Stopping is symmetric
critical_means_boundary <- function(ns, mu, prior_var = 1, bf_crit = 3, step = 1, symmetric = TRUE) {
  # Define the range for the normal distribution
  y_axis <- seq(-5, 5, 0.01)
  x_norm <- dnorm(y_axis, mean = mu, sd = sqrt(prior_var) / sqrt(step))
  if (symmetric) {
    fn <- "figures/sym-critical-means.pdf"
  } else {
    fn <- "figures/asym-critical-means.pdf"
  }
  pdf(fn)
  plot(0, 0,
    xlim = c(-1, 50), ylim = c(-2, 2),
    main = bquote("Decision Boundaries for " * H[0] * " and " * H[1]),
    ylab = bquote("Estimated Mean " * bar(y)), xlab = "Sample Size",
    type = "n"
  )
  y_crit_h1 <- get_y_by_bf(ns, prior_var, bf_crit)
  y_crit_h0 <- get_y_by_bf(ns, prior_var, 1 / bf_crit)
  # plot mu with sem
  # abline(h = mu, col = "black", lwd = 2)

  lines(ns, mu + 1.96 * sqrt(prior_var) / sqrt(ns),
    col = "orange", lty = 2, lwd = 2
  )
  lines(ns, mu - 1.96 * sqrt(prior_var) / sqrt(ns),
    col = "orange", lty = 2, lwd = 2
  )
  # Normal distribution around mu
  abline(v = 0, lwd = 2, col = "black")
  lines(-x_norm, y_axis, col = "black", lwd = 2)

  # Find the corresponding y-values where n = step
  # and plot the segment from the normal distribution to the decision boundary
  y1_crit_h0 <- get_y_by_bf(step, prior_var, 1 / bf_crit)
  y2_crit_h0 <- -get_y_by_bf(step, prior_var, 1 / bf_crit)
  segments(step, y1_crit_h0, -dnorm(y1_crit_h0, mu, sqrt(prior_var) / sqrt(step)), y1_crit_h0, col = "black", lty = 2)
  segments(step, y2_crit_h0, -dnorm(y2_crit_h0, mu, sqrt(prior_var) / sqrt(step)), y2_crit_h0, col = "black", lty = 2)
  if (symmetric) {
    y1_crit_h1 <- get_y_by_bf(step, prior_var, bf_crit)
    y2_crit_h1 <- -get_y_by_bf(step, prior_var, bf_crit)
    segments(step, y1_crit_h1, -dnorm(y1_crit_h1, mu, sqrt(prior_var) / sqrt(step)), y1_crit_h1, col = "black", lty = 2)
    segments(step, y2_crit_h1, -dnorm(y2_crit_h1, mu, sqrt(prior_var) / sqrt(step)), y2_crit_h1, col = "black", lty = 2)
  }
  # plot h0 decision boundaries
  lines(ns, y_crit_h0, col = "#59B3E6", lwd = 2)
  lines(ns, -y_crit_h0, col = "#59B3E6", lwd = 2)
  if (symmetric) {
    # plot h1 decision boundaries
    lines(ns, y_crit_h1, col = "#CD1076", lwd = 2)
    lines(ns, -y_crit_h1, col = "#CD1076", lwd = 2)
  }

  legend("topright",
    legend = c(
      "Decision for H0", if (symmetric) {
        "Decision for H1"
      },
      bquote("Normal Distribution around " * mu == .(mu)),
      "95% Confidence Interval",
      "Random Walk"
    ),
    fill = c("#59B3E6", if (symmetric) {
      "#CD1076"
    }, "black", "orange", "#009980")
  )
  # add random walk by comparing mean of normal distributed mean with critical values
  # for decision making
  if (symmetric) {
    y1_crit_h0 <- c()
    y2_crit_h0 <- c()
    y1_crit_h1 <- c()
    y2_crit_h1 <- c()
    for (n in ns) {
      y1_crit_h0 <- c(y1_crit_h0, get_y_by_bf(n, prior_var, 1 / bf_crit))
      y2_crit_h0 <- c(y2_crit_h0, -get_y_by_bf(n, prior_var, 1 / bf_crit))
      y1_crit_h1 <- c(y1_crit_h1, get_y_by_bf(n, prior_var, bf_crit))
      y2_crit_h1 <- c(y2_crit_h1, -get_y_by_bf(n, prior_var, bf_crit))
    }
    nan_to_inf <- function(x, sign = "+") {
      if (is.nan(x) && sign == "+") {
        Inf
      } else if (is.nan(x) && sign == "-") {
        -Inf
      } else {
        x
      }
    }
    y1_crit_h0 <- sapply(y1_crit_h0, function(x) {
      nan_to_inf(x, sign = "-")
    })
    y2_crit_h0 <- sapply(y2_crit_h0, function(x) {
      nan_to_inf(x, sign = "+")
    })
    y1_crit_h1 <- sapply(y1_crit_h1, function(x) {
      nan_to_inf(x, sign = "+")
    })
    y2_crit_h1 <- sapply(y2_crit_h1, function(x) {
      nan_to_inf(x, sign = "-")
    })
    # plot random walk
    n <- 1
    y <- rnorm(1, mu, sqrt(prior_var) / sqrt(n))
    y_mean <- mean(y)
    y_list <- c(y)
    mean_y_list <- c(y_mean)
    while (!(((y1_crit_h0[n] != y2_crit_h0[n]) && (y_mean <= y1_crit_h0[n] && y_mean >= y2_crit_h0[n])) || (y1_crit_h0[n] == y2_crit_h0[n] && y2_crit_h0[n] == y_mean) || (y_mean >= y1_crit_h1[n] || y_mean <= y2_crit_h1[n]))) {
      if (n < length(y1_crit_h0)) {
        n <- n + 1
      } else {
        print(paste("Aborted random walk after", n, "steps."))
        break
      }
      y <- rnorm(1, mu, sqrt(prior_var) / sqrt(n))
      y_list <- c(y_list, y)
      y_mean <- mean(y_list)
      mean_y_list <- c(mean_y_list, y_mean)
    }
    print(mean_y_list[n])
    print(n)
    print(y1_crit_h1[n])
    print(y2_crit_h1[n])
    points(n, y_mean, pch = 4, col = "#009980", cex = 2)
    lines(seq_along(mean_y_list), mean_y_list, col = "#009980", lwd = 2)
  } else {
    y1_crit_h0 <- c()
    y2_crit_h0 <- c()
    for (n in ns) {
      y1_crit_h0 <- c(y1_crit_h0, get_y_by_bf(n, prior_var, 1 / bf_crit))
      y2_crit_h0 <- c(y2_crit_h0, -get_y_by_bf(n, prior_var, 1 / bf_crit))
    }
    nan_to_inf <- function(x, sign = "+") {
      if (is.nan(x) && sign == "+") {
        Inf
      } else if (is.nan(x) && sign == "-") {
        -Inf
      } else {
        x
      }
    }
    y1_crit_h0 <- sapply(y1_crit_h0, function(x) {
      nan_to_inf(x, sign = "-")
    })
    y2_crit_h0 <- sapply(y2_crit_h0, function(x) {
      nan_to_inf(x, sign = "+")
    })
    # plot random walk
    n <- 1
    y <- rnorm(1, mu, sqrt(prior_var) / sqrt(n))
    y_mean <- mean(y)
    y_list <- c(y)
    mean_y_list <- c(y_mean)
    while (!(((y1_crit_h0[n] != y2_crit_h0[n]) && (y_mean <= y1_crit_h0[n] && y_mean >= y2_crit_h0[n])) || (y1_crit_h0[n] == y2_crit_h0[n] && y2_crit_h0[n] == y_mean))) {
      if (n < length(y1_crit_h0)) {
        n <- n + 1
      } else {
        print(paste("Aborted random walk after", n, "steps."))
        break
      }
      y <- rnorm(1, mu, sqrt(prior_var) / sqrt(n))
      y_list <- c(y_list, y)
      y_mean <- mean(y_list)
      mean_y_list <- c(mean_y_list, y_mean)
    }
    print(mean_y_list[n])
    print(n)
    points(n, y_mean, pch = 4, col = "#009980", cex = 2)
    lines(seq_along(mean_y_list), mean_y_list, col = "#009980", lwd = 2)
  }
  dev.off()
}

#################################
### PLOT DECISION PROBABILITY ###
#################################
# plot the decision probability for H0 given mu and n_start
# case: "asymmetrical", "symmetrical", "both"
# n_start: number of starts for the random walk
# approx: use approximation, default is FALSE
decision_prob_curve <- function(case = "both", n_start, approx = FALSE) {
  # plot the decision probability for H0 given mu
  pdf(paste("figures/", case, if (approx) {
    "-approx"
  }, if (n_start > 1) {
    paste("-n_start-", n_start, sep = "")
  }, "-decision-probability.pdf", sep = ""))
  main_text <- if (case == "both") {
    " rules"
  } else {
    " rule"
  }
  approx_text <- if (approx) {
    " with approximation"
  } else {
    ""
  }
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = bquote("Decision probability for " * .(case) * .(main_text) * .(approx_text)),
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

  # approximation
  if (approx) {
    # asymmetrical
    if (case == "asymmetrical") {
      asym_h0_probs <- c()
      is_ph0_smaller_50 <- FALSE
      for (mu in big_sim_mus) {
        asym_h0_prob <- tail(asym_decision_probability(mu, 1, 3, n_start), n = 1)
        asym_h0_probs <- c(
          asym_h0_probs,
          asym_h0_prob
        )
        if (!is_ph0_smaller_50 && asym_h0_prob < 0.5) {
          is_ph0_smaller_50 <- TRUE
          print(paste("P(H0) > 0.5 for mu =", mu))
        }
      }
      lines(big_sim_mus, asym_h0_probs, col = "#009980", lwd = 2)
      legend("topright",
        legend = c("Approximation", "Simulation"),
        fill = c("#009980", "orange")
      )
    } else if (case == "symmetrical") {
      sym_h0_probs <- c()
      is_ph0_smaller_50 <- FALSE
      for (mu in big_sim_mus) {
        sym_h0_prob <- tail(sym_decision_probability(mu, 1, 3, n_start), n = 1)
        sym_h0_probs <- c(
          sym_h0_probs,
          sym_h0_prob
        )
        if (!is_ph0_smaller_50 && sym_h0_prob < 0.5) {
          is_ph0_smaller_50 <- TRUE
          print(paste("P(H0) > 0.5 for mu =", mu))
        }
      }
      lines(big_sim_mus, sym_h0_probs, col = "blue", lwd = 2)
      legend("topright",
        legend = c("Approximation", "Simulation"),
        fill = c("blue", "#FF99FF")
      )
    } else {
      print("Approximation is not implemented for both cases. Its not a bug, its a feature.")
    }
  } else {
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
  }
  dev.off()
}
# Setup to save the plot as pdf
sim_diff_n_start <- function() {
  # import library
  library("RColorBrewer")
  # config
  colors <- heat.colors(15)
  pdf("figures/dbh0-diff-all-n-starts.pdf")
  # Increase right margin to make space for the color bar
  par(mar = c(4, 6, 4, 6), xpd = TRUE) # Increase right margin and allow plotting outside the plot

  # Create the plot
  plot(0, 0,
    type = "l",
    lwd = 2, xlab = bquote(mu), ylab = bquote(P(H[0])), main = bquote(P(H[0]) * " for difference simulation " - " approximation"),
    xlim = c(0, 1), ylim = c(-0.5, 0.5)
  )


  labels <- 1:15

  for (n_start in 1:15) {
    con <- dbConnect(duckdb(), "data/hacking-bayes.duckdb")
    sym <- data.frame()
    sym <- dbGetQuery(
      con, "SELECT *
                            FROM bf_decision_threshold
                            WHERE symmetrical = TRUE
                            AND trial_start = ?",
      list(n_start)
    )
    dbDisconnect(con)
    h0_data <- sym[sym[, 1] == 0, , ]
    h1_data <- sym[sym[, 1] == 1, , ]

    h0_prob <- c()
    for (mu in big_sim_mus) {
      h0_count <- length(h0_data[h0_data[, 3] == mu, 2])
      h0_prob <- c(h0_prob, h0_count / repetitions)
    }

    sym_h0_prob <- c()
    for (mu in big_sim_mus) {
      sym_h0_prob <- c(sym_h0_prob, tail(sym_decision_probability(mu, 1, 3, n_start), n = 1))
    }

    diff <- h0_prob - sym_h0_prob # Difference
    lines(big_sim_mus, diff, col = colors[n_start], lwd = 2)
  }

  par(xpd = FALSE)
  abline(h = 0, col = "black", lty = 2)
  # Enablin plotting outside the plot area
  par(xpd = TRUE)
  # Add the discrete color cartridge
  n_colors <- length(colors)
  y_positions <- seq(-0.5, 0.5, length.out = n_colors + 1) # Positions for discrete color blocks

  # Plot the discrete color blocks
  for (i in 1:n_colors) {
    rect(1.08, y_positions[i], 1.13, y_positions[i + 1], col = colors[i], border = NA)
  }

  # Add text labels next to the color blocks
  for (i in 1:n_colors) {
    text(1.15, (y_positions[i] + y_positions[i + 1]) / 2, labels[i], adj = 0) # Adjust text position
  }
  # Save and close the pdf device
  dev.off()
}
