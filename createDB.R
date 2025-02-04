# clear workspace
rm(list = ls())
setwd("/home/linus/git/hacking-bayes")
# load libraries
library("duckdb")

# create the database
# con <- dbConnect(duckdb(), dbdir = "data/results.duckdb", read_only = FALSE)
# dbDisconnect(con)

# migrate the tables from realistic_sym_simulation.db to results.duckdb
con <- dbConnect(duckdb(), dbdir = "data/results.duckdb", read_only = FALSE)
dbExecute(con, "ATTACH DATABASE 'data/realistic_sym_simulation.db' AS sym_db")
dbExecute(con, "CREATE TABLE cauchy_sym AS SELECT * FROM sym_db.realistic_sym_simulation")
dbExecute(con, "CREATE TABLE cauchy_prior AS SELECT * FROM sym_db.cauchy_fun")
dbDisconnect(con)

# OLD CODE
# migrate RDS files to results.duckdb all to the table bf_decision_thresholds
con <- dbConnect(duckdb(), dbdir = "data/results.duckdb", read_only = FALSE)
dbExecute(con, "DROP TABLE IF EXISTS bf_decision_threshold")
dbExecute(con, "CREATE TABLE bf_decision_threshold (decision INTEGER, stop_count INTEGER, mu DOUBLE, trial_start INTEGER, symmetrical BOOLEAN)")
for (i in c(1:15, 30)) {
    # Load RDS
    rdsfile <- readRDS(paste0("data/from_", i, "_asym_simulation_data", sep = ""))
    # Add column for symmetrical
    rdsfile$trial_start <- i
    rdsfile$symmetrical <- FALSE
    setnames(rdsfile, c("decision", "stop_count", "mu", "trial_start", "symmetrical"))
    # Write to duckdb
    dbWriteTable(con, "bf_decision_threshold", rdsfile, append = TRUE)
}
for (i in c(1:15, 30)) {
    # Load RDS
    rdsfile <- readRDS(paste0("data/from_", i, "_sym_simulation_data", sep = ""))
    # Add column for symmetrical
    rdsfile$trial_start <- i
    rdsfile$Symmetrical <- TRUE
    setnames(rdsfile, c("decision", "stop_count", "mu", "trial_start", "symmetrical"))
    # Write to duckdb
    dbWriteTable(con, "bf_decision_threshold", rdsfile, append = TRUE)
}
dbDisconnect(con)
