setwd("/home/linus/git/hacking-bayes")
# clear workspace
rm(list = ls())
# Load required libraries
library(BayesFactor)
library(foreach)
library(doParallel)

# Register parallel backend
registerDoParallel(cores = detectCores())

# Define the parallelized simulation function
cauchy_plot_parallel <- function(repetitions, max_n, mus, r) {
    begin_time <- Sys.time()
    total_rows <- repetitions * length(mus) * (max_n - 2) # Total number of rows
    temp_files <- vector("list", length(mus) * repetitions) # List to store temporary file names

    # Function to perform the simulation
    simulate_for_mu_repetition <- function(mu_value, rep_id, temp_file) {
        print(paste("Simulation starts for mu =", mu_value, ", repetition =", rep_id))
        start <- Sys.time()
        data <- rnorm(2, mean = mu_value, sd = 1) # Initial data
        bf_list <- numeric(max_n - 2) # To store Bayes Factors
        n_list <- integer(max_n - 2) # To store sample sizes
        row_id <- 1 # Track the row id in the file

        # Collect the Bayes Factor at every step
        for (n in 3:max_n) { # Start from 3 since we already have 2 data points
            data <- c(data, rnorm(1, mean = mu_value, sd = 1))
            bf <- extractBF(ttestBF(data, mu = 0, rscale = r))$bf
            bf_list[row_id] <- bf
            n_list[row_id] <- n
            row_id <- row_id + 1
        }

        # Write results to the temporary file
        write.table(data.frame(repetition = rep_id, mu = mu_value, n = n_list, bf = bf_list),
            file = temp_file, row.names = FALSE, col.names = !file.exists(temp_file), append = TRUE, sep = ","
        )
        end <- Sys.time()
        print(paste("Simulation done for mu =", mu_value, ", repetition =", rep_id, "in", end - start))
    }

    # Parallelize the simulation
    foreach(mu = mus, .combine = "c", .packages = c("BayesFactor")) %:%
        foreach(rep_id = 1:repetitions, .combine = "c") %dopar% {
            temp_file <- tempfile(fileext = ".csv")
            simulate_for_mu_repetition(mu, rep_id, temp_file)
            temp_file
        } -> temp_files

    # Combine all temporary files into one CSV file
    final_file <- "data/cauchy_prior_simulation.csv"
    file_conn <- file(final_file, open = "wt")
    write.csv(data.frame(repetition = integer(), mu = double(), n = integer(), bf = double()),
        file = file_conn, row.names = FALSE
    )

    for (temp_file in temp_files) {
        file_content <- readLines(temp_file)
        writeLines(file_content[-1], file_conn) # Skip header for all files except the first
        unlink(temp_file) # Remove the temporary file
    }
    close(file_conn)
    end_time <- Sys.time()
    print(paste("Everything done in", end_time - begin_time))
}

# Example usage
repetitions <- 20000
max_n <- 100
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
r <- 1 / sqrt(2)

# Run the parallelized simulation
cauchy_plot_parallel(repetitions, max_n, mus, r)
print(paste("Total time taken:", end_time - begin_time))

# Stop the parallel backend
stopImplicitCluster()
