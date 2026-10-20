# Simulating the Bayesian Optional Stopping procedure 
# with different cauchy priors and storing each sampled 
# data point and calculated Bayes Factor.
setwd(".")
# Clear workspace
rm(list = ls())
# Load required libraries
library(BayesFactor)
library(foreach)
library(doParallel)
library(data.table)
library(duckdb)

# Register parallel backend
registerDoParallel(cores = detectCores() - 1)

#'@description Calculate one Bayesian Optional Stopping procedure 
#' and store results until maximum sample size is reached.
#' 
#' @param mu_value true effect size
#' @param repl_id id of the current replication
#' @param max_n maximum number of samples to collect
#' @return data.table with sample size and corresponding Bayes Factor
simulate_for_mu_repetition <- function(mu_value, repl_id, max_n) {
    print(paste("Simulation starts for mu =", mu_value, ", repetition =", repl_id))
    start <- Sys.time()
    data <- rnorm(2, mean = mu_value, sd = 1) # Initial data
    results <- data.table(repetition = numeric(max_n - 2), mu = numeric(max_n - 2), n = numeric(max_n - 2), bf = numeric(max_n - 2))
    id <- 1 # Track the row id in the file

    # Collect the Bayes Factor at every step
    # (loop variables must not be called n or bf: inside `:=` those names refer to the columns)
    for (sample_size in 3:max_n) { # Start from 3 since we already have 2 data points
        data <- c(data, rnorm(1, mean = mu_value, sd = 1))
        bf_value <- extractBF(ttestBF(data, mu = 0, rscale = r))$bf
        results[id, `:=`(repetition = repl_id, mu = mu_value, n = sample_size, bf = bf_value)]
        id <- id + 1
    }
    end <- Sys.time()
    # print(results)
    print(paste("Simulation done for mu =", mu_value, ", repetition =", repl_id, "in", end - start))
    return(results)
}

##########################
####### SIMULATION #######
##########################
repetitions <- 20000
# Adapt chunk_size depending on your memory capacity for parallelization
chunk_size <- 2000
chunks <- 1:(repetitions / chunk_size)
max_n <- 100
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
r <- 1 / sqrt(2)

begin_time <- Sys.time()
total_rows <- repetitions * length(mus) * (max_n - 2) # Total number of rows
temp_files <- vector("list", length(mus) * repetitions) # List to store temporary file names

# Open database connection and create corresponding table
con <- dbConnect(duckdb(), dbdir = "data/hacking-bayes.duckdb")
dbExecute(con, "DROP TABLE IF EXISTS cauchy_prior")
dbExecute(con, "CREATE TABLE cauchy_prior (repetition DOUBLE, mu DOUBLE, n DOUBLE, bf DOUBLE)")

# Parallelize for each true effect size and replication and specified chunk size
foreach(mu = mus) %do% {
    for (chunk in chunks) {
        chunk_start <- 1 + (chunk - 1) * chunk_size
        chunk_end <- chunk * chunk_size
        print(paste("Starting chunk", chunk, "from", chunk_start, "to", chunk_end))
        rbindlist(foreach(rep_id = chunk_start:chunk_end) %dopar% {
            simulate_for_mu_repetition(mu, rep_id, max_n)
        }) -> results
        # Write results in the database table
        dbWriteTable(con, "cauchy_prior", results, append = TRUE)
    }
}
dbDisconnect(con)
# Stop the parallel backend
stopImplicitCluster()

end_time <- Sys.time()
# Run the parallelized simulation
print(paste("Total time taken:", end_time - begin_time))
