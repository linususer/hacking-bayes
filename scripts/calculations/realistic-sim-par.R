# Set working directory
setwd("/home/linus/git/hacking-bayes")

# Clear workspace
rm(list = ls())
gc()

# Import necessary libraries
library(BayesFactor)
library(foreach)
library(doParallel)

# Register parallel backend
registerDoParallel(cores = detectCores())

# Function to perform the simulation
simulate_bayes_factor <- function(mu, r, BF_crit, repetitions, temp_file) {
  start <- Sys.time()
  # Open a file connection
  file_conn <- file(temp_file, open = "wt")

  # Write the column headers
  write.csv(data.frame(Decision_H0_H1 = integer(), Stop_Count = integer(), Mu = double(), BF_crit = double(), r = double()),
    file = file_conn, row.names = FALSE
  )

  # Simulation
  for (i in 1:repetitions) {
    x <- rnorm(2, mean = mu, sd = 1)
    BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
    stop_count <- 1
    # asym
    while (is.na(BF) || is.null(BF) || (BF > (1 /BF_crit) && stop_count < 250)) { # sym: && BF < BF_crit)) {
      x <- c(x, rnorm(1, mean = mu, sd = 1))
      BF <- extractBF(ttestBF(x, mu = 0, r = r))$bf
      stop_count <- stop_count + 1
    }
    decision <- ifelse(BF < BF_crit, 0, 1)
    result <- data.frame(Decision_H0_H1 = decision, Stop_Count = stop_count, Mu = mu, BF_crit = BF_crit, r = r)

    # Write result to the file
    write.table(result, file = file_conn, row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  }
  close(file_conn)
  end <- Sys.time()
  print(paste("Simulation done for mu =", mu, "BF_crit =", BF_crit, "r =", r, "in", end - start))
}

# Initialize variables
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3,6,10)
repetitions <- 20

# Perform parallel simulations
temp_files <- foreach(BF_crit = BF_crits) %:%
  foreach(r = r_vals) %:%
  foreach(mu = big_sim_mus) %dopar% {
    temp_file <- paste0(tempfile(), ".csv")
    print(paste("Starting simulation for mu =", mu, "BF_crit =", BF_crit, "r =", r))
    simulate_bayes_factor(mu, r, BF_crit, repetitions, temp_file)
    return(temp_file)
  }

final_file <- "data/realistic_asym_simulation.csv"
file_conn <- file(final_file, open = "wt")
write.csv(data.frame(Decision_H0_H1 = integer(), Stop_Count = integer(), Mu = double(), BF_crit = double(), r = double()),
  file = file_conn, row.names = FALSE
)

temp_files <- unlist(temp_files)
for (temp_file in temp_files) {
  file_content <- readLines(temp_file)
  writeLines(file_content[-1], file_conn)
  unlink(temp_file)
}
close(file_conn)

# Stop the parallel backend
stopImplicitCluster()
print("DONE")
