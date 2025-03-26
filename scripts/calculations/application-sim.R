# A simluation to create an application example how the misspecification 
# of each parameter can be used to tweak results in Bayesian Optional Stopping
# even though the data is drawn from the same distribution.
setwd(".")
# Clear workspace
rm(list = ls())
# Load libraries
library("data.table")
library("duckdb")
# Import Bayes Factor functions
source("bayes-factor-functions.R")

#' @description Calculate an applicable example of misspecified params 
#' with Bayesian Optional Stopping drawn from a normal distribution.
#'
#' @param mu true effect size
#' @param fixed_sigma standard deviation of the specified prior
#' @param prior_var true, known variance of the drawn data
#' @param bf_crit critical Bayes Factor threshold for H1
#' @param n_start sample to start with Optional Stopping
#' @param n_end sample to stop with Optional Stopping
#' @return resulting data.table of 20000 Optional Stopping decisions
#' @note for reference see Szillat (2024), p. 58ff. 
#' @references Szillat, Linus R.P. (2024). Hacking Bayes Factors with Optional Stopping 
#' [Unpublished bachelor's thesis]. University of Tübingen. Also can be found at linus-szillat.de/thesis.
applicationExampleSim <- function(mu, fixed_sigma, prior_var, bf_crit, n_start, n_end) {
    results <- vector("list", length = 20000)
    repetitions <- 20000
    # reset h0/h1_counts
    h0_count <- 0
    h1_count <- 0
    # run the simulation 20000 times
    for (i in 1:repetitions) {
        h1 <- rnorm(n_start, mean = mu, sd = fixed_sigma)
        bf <- bf01(bf10_general(length(h1), mean(h1), var(h1), prior_var))
        if (is.na(bf)) {
            bf <- 1.33
        }
        # run the simulation until the bayes factor is between 1/3 and 3
        stop_count <- 1
        while (bf < bf_crit && bf > 1 / bf_crit && stop_count < n_end) {
            h1 <- c(h1, rnorm(1, mean = mu, sd = fixed_sigma))
            bf <- bf01(bf10_general(length(h1), mean(h1), var(h1), prior_var))
            stop_count <- stop_count + 1
        }
        # count the number of times the null hypothesis is true
        # and save the results for plotting
        if (bf > bf_crit) {
            h0_count <- h0_count + 1
            results[[h0_count + h1_count]] <- list(0, stop_count, mu, n_start, n_end)
        } else {
            h1_count <- h1_count + 1
            results[[h0_count + h1_count]] <- list(1, stop_count, mu, n_start, n_end)
        }
        # print the results
        # print(paste("Decision: ", bf, "Stop count: ", stop_count))
    }
    # print the results
    print(paste("Symmetrical simulation is done for mu =", mu))
    print(paste("H0: ", h0_count))
    print(paste("H1: ", h1_count))
    # Export results as data frame
    results <- rbindlist(results)
    setnames(results, c("decision", "stop_count", "mu", "trial_start", "trial_end"), )
    results
}


# Application Example to maximize P('H0')
# Set configuration

mu <- 0.1
fixed_sigma <- 2
prior_var <- 3^2
bf_crit <- 3
ns <- get_n_by_bf(0.1, fixed_sigma^2, prior_var, bf_crit)
n_start <- if (is.nan(ns[1])) 1 else round(ns[2])
n_end <- if (is.nan(ns[2])) Inf else round(ns[1])

# Calculate results and store them in a duckDB database
results <- applicationExampleSim(mu, fixed_sigma, prior_var, bf_crit, n_start, n_end)
con <- dbConnect(duckdb(), dbdir = "data/results.duckdb", read_only = FALSE)
dbExecute(con, "DROP TABLE IF EXISTS application_example")
dbWriteTable(con, "application_example", results)
dbDisconnect(con)

# Application example to maximize P('H1')
# Set configuration
mu <- 0.1
fixed_sigma <- 1
prior_var <- 0.15^2
bf_crit <- 3
ns <- get_n_by_bf(0.1, fixed_sigma^2, prior_var, bf_crit)
n_start <- if (is.nan(ns[1])) 1 else round(ns[2])
n_end <- if (is.nan(ns[2])) Inf else round(ns[1])

# Calculate results and store them in a duckDB database
results <- applicationExampleSim(mu, fixed_sigma, prior_var, bf_crit, n_start, n_end)
con <- dbConnect(duckdb(), dbdir = "data/results.duckdb", read_only = FALSE)
dbAppendTable(con, "application_example", results)
dbDisconnect(con)
