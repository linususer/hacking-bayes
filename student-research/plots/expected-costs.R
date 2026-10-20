# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(duckdb)
source("shared/define_colors.R")
source("shared/queries.R")
source("student-research/plots/expected-costs-summaries.R")
source("student-research/plots/expected-costs-plots.R")
# Initialise params
big_sim_mus <- seq(0, 1, 0.01)
r_vals <- c(0.5, 1, 2) / sqrt(2)
BF_crits <- c(3, 6, 10)
repetitions <- 20000

db_file <- "data/hacking-bayes.duckdb"
db_large_effects_file <- "data/large_effects.duckdb"


#####################################
####### Expected cost summaries #########
#####################################


#########################
######## CONFIG #########
#########################

# Cost matrix 1: symmetrical costs and indecisive outcomes cost > Type I & Type II error cost
cost_matrix1 <- matrix(
    c(
        0, 0.5, 1,
        0.5, 0, 1
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)
# Cost matrix 2: symmetrical costs and indecisive outcomes = Type 1 & Type II error cost
cost_matrix2 <- matrix(
    c(
        0, 1, 1,
        1, 0, 1
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)
# Cost matrix 3: symmetrical costs and indecisive outcomes cost < Type I & Type II error cost
cost_matrix3 <- matrix(
    c(
        0, 1, 0.5,
        1, 0, 0.5
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)

# Plots for expected costs in general
plot_expected_costs(3, r_vals[2], scenario_name = "default")
plot_expected_costs(3, r_vals[2], cost_matrix1, scenario_name = "sym-ind-ge-error")
plot_expected_costs(3, r_vals[2], cost_matrix2, scenario_name = "sym-ind-e-error")
plot_expected_costs(3, r_vals[2], cost_matrix3, scenario_name = "sym-ind-le-error")
# Plots for BF_crit
plot_cost_by_bf(scenario_name = "default")
plot_cost_by_bf(cost_matrix1, scenario_name = "sym-ind-ge-error")
plot_cost_by_bf(cost_matrix2, scenario_name = "sym-ind-e-error")
plot_cost_by_bf(cost_matrix3, scenario_name = "sym-ind-le-error")


# Cost matrix 4:
cost_matrix4 <- matrix(
    c(
        0, 1, 0.25,
        1, 0, 0.25
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)

plot_expected_costs(3, r_vals[2], cost_matrix4, scenario_name = "sym-ind-bfcrit3")
plot_cost_by_bf(cost_matrix4, scenario_name = "sym-ind-bfcrit3")

cost_matrix5 <- matrix(
    c(
        0, 1, 1 / 7,
        1, 0, 1 / 7
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)

cost_matrix6 <- matrix(
    c(
        0, 1, 1 / 11,
        1, 0, 1 / 11
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("H0_true", "H1_true"),
        c("dec_H0", "dec_H1", "dec_ind")
    )
)

plot_expected_costs(6, r_vals[2], cost_matrix5, scenario_name = "sym-ind-bfcrit6")
plot_cost_by_bf(cost_matrix5, scenario_name = "sym-ind-bfcrit6")
plot_expected_costs(10, r_vals[2], cost_matrix6, scenario_name = "sym-ind-bfcrit10")
plot_cost_by_bf(cost_matrix6, scenario_name = "sym-ind-bfcrit10")


###########################################
############ WEIGHT COSTS #################
###########################################


# NEW: map each BF_crit to its purpose-built cost matrix
bf_crit_cost_matrices <- list(
    `3`  = cost_matrix4,
    `6`  = cost_matrix5,
    `10` = cost_matrix6
)

con <- dbConnect(duckdb(), db_file)
for (bfc in c(3)) { # , 6, 10)) {
    cm <- bf_crit_cost_matrices[[as.character(bfc)]]
    for (dtype in c("opt_stop", "stop_avg", "same_dist")) {
        for (mode in c("interval", "threshold")) { # NEW: both extent modes, every call
            plot_prior_distr_cost_summary(prior_distr_cost_summary(con,
                bf_crit = bfc, r_val = r_vals[2],
                cost_matrix = cm, df_type = dtype,
                grid_extent_mode = mode,
                mass_interval = c(-1, 1),
                tail_threshold = 10e-5,
                scenario_name = sprintf("sym-ind-bfcrit%s", bfc),
                summary_path = sprintf(
                    "student-research/summaries/prior-weighted-cost_%s-sym-ind-bfcrit%s-bf%s-%s.txt",
                    dtype, bfc, bfc, mode
                )
            ))
        }
    }
}
dbDisconnect(con)

# stop avg
# 37          1.80   0.08142341
# 38          1.85   0.07800663
# 39          1.90   0.07458985
# 40          1.95   0.07117308
# 41          2.00   0.06775630
# 42          2.05   0.06433952
# 43          2.10   0.06092274
# 44          2.15   0.05750597
# 45          2.20   0.05408919
# 46          2.25   0.05067241


con <- dbConnect(duckdb(), db_file)
con_large <- dbConnect(duckdb(), db_large_effects_file)

all_summaries   <- list() # one summary for all results
per_bf_summaries <- list() # one summary per bf crit

for (bfc in c(3, 6, 10)) {
  cm <- bf_crit_cost_matrices[[as.character(bfc)]]
  bf_key <- as.character(bfc)
  per_bf_summaries[[bf_key]] <- list()

  for (dtype in c("opt_stop", "stop_avg", "same_dist")) {
    for (bound in c("upper", "lower")) {

      summary <- summary_expected_costs(
        con,
        bf_crit = bfc, r_val = r_vals[2],
        bound = bound,
        cost_matrix = cm, df_type = dtype,
        mass_interval = c(-3, 3),
        con_large = con_large
      )
      summary_df <- as.data.frame(summary)

      all_summaries[[length(all_summaries) + 1]] <- summary_df
      per_bf_summaries[[bf_key]][[length(per_bf_summaries[[bf_key]]) + 1]] <- summary_df
    }
  }
  df_bf <- do.call(rbind, per_bf_summaries[[bf_key]])
  write.table(
    df_bf,
    sprintf("student-research/summaries/expected_costs_bf-%s.txt", bfc),
    sep = "\t", row.names = FALSE, col.names = TRUE
  )
}

df_all <- do.call(rbind, all_summaries)
write.table(
  df_all,
  "student-research/summaries/expected_costs_all.txt",
  sep = "\t", row.names = FALSE, col.names = TRUE
)
dbDisconnect(con)
dbDisconnect(con_large, shutdown = TRUE)

# Plot the worst-case cost per bin from the saved summaries
for (bfc in c(3, 6, 10)) {
    summary_title <- sprintf("student-research/summaries/expected_costs_bf-%s.txt", bfc)
    for (dtype in c("opt_stop", "stop_avg", "same_dist")) {
        for (bound in c("upper", "lower")) {
            plot_cost_bound(summary_title, df_type = dtype, bound = bound)
        }
    }
}

summary_all <- read.table(
  "student-research/summaries/expected_costs_all.txt",
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, quote = "\""
)

## overview plot for everything


plot_overview_expected_costs(summary_all)
