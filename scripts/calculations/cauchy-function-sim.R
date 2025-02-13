setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# Load required libraries
library(BayesFactor)
library(foreach)
library(doParallel)
library(data.table)
library(duckdb)

# Register parallel backend
registerDoParallel(cores = detectCores() - 1)

# Define the parallelized simulation function
# Function to perform the simulation
simulate_for_mu_repetition <- function(mu_value, repl_id) {
    print(paste("Simulation starts for mu =", mu_value, ", repetition =", repl_id))
    start <- Sys.time()
    data <- rnorm(2, mean = mu_value, sd = 1) # Initial data
    results <- data.table(rep_id = integer(max_n - 2), mu = numeric(max_n - 2), sample_size = integer(max_n - 2), bayes_factor = numeric(max_n - 2))
    id <- 1 # Track the row id in the file

    # Collect the Bayes Factor at every step
    for (n in 3:max_n) { # Start from 3 since we already have 2 data points
        data <- c(data, rnorm(1, mean = mu_value, sd = 1))
        bf <- extractBF(ttestBF(data, mu = 0, rscale = r))$bf
        results[id, `:=`(rep_id = repl_id, mu = mu_value, sample_size = n, bayes_factor = bf)]
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
# Example usage
repetitions <- 20000
# Adapt chunk_size depending on your memory capacity
chunk_size <- 2000
chunks <- 1:(repetitions / chunk_size)
max_n <- 100
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
r <- 1 / sqrt(2)

begin_time <- Sys.time()
total_rows <- repetitions * length(mus) * (max_n - 2) # Total number of rows
temp_files <- vector("list", length(mus) * repetitions) # List to store temporary file names

# Parallelize the simulation
con <- dbConnect(duckdb(), dbdir = "data/results.duckdb")
dbExecute(con, "DROP TABLE IF EXISTS cauchy_prior_simulation")
dbExecute(con, "CREATE TABLE cauchy_prior_simulation (rep_id INTEGER, mu DOUBLE, sample_size INTEGER, bayes_factor DOUBLE)")

foreach(mu = mus) %do% {
    for (chunk in chunks) {
        chunk_start <- 1 + (chunk - 1) * chunk_size
        chunk_end <- chunk * chunk_size
        print(paste("Starting chunk", chunk, "from", chunk_start, "to", chunk_end))
        rbindlist(foreach(rep_id = chunk_start:chunk_end) %dopar% {
            simulate_for_mu_repetition(mu, rep_id)
        }) -> results
        dbWriteTable(con, "cauchy_prior_simulation", results, append = TRUE)
    }
}
dbDisconnect(con)
# Stop the parallel backend
stopImplicitCluster()

end_time <- Sys.time()
# Run the parallelized simulation
print(paste("Total time taken:", end_time - begin_time))
