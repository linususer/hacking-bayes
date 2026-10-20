# Rebuilds the figures of all projects, optionally after rebuilding their summaries.
#
#   Rscript run-all.R              figures only, from the committed summaries (no database needed)
#   Rscript run-all.R summarise    first the summaries from data/*.duckdb, then the figures
#
# The simulations in */calculations/ are never run by this script.
# Every script runs in its own R process because each one clears its workspace.

summarise_scripts <- c(
  "thesis/summarise/realistic-sim.R",
  "thesis/summarise/decision-prob-bf-crit.R",
  "thesis/summarise/optional-stopping.R",
  "lab-project/summarise/fixed-size.R",
  "student-research/summarise/expected-costs.R"
)

plot_scripts <- c(
  "thesis/plots/binom-example.R",
  "thesis/plots/catch-up-effect-plot-functions.R",
  "thesis/plots/realistic-sim.R",
  "thesis/plots/decision-prob-bf-crit.R",
  "thesis/plots/optional-stopping.R",
  "lab-project/plots/fixed-size.R",
  "lab-project/plots/random_walk_priors.R",
  "student-research/plots/expected-costs.R",
  "student-research/plots/e-values.R"
)

run_script <- function(script) {
  cat("==>", script, "\n")
  status <- system2("Rscript", script)
  if (status != 0) stop(script, " failed with exit status ", status)
}

if ("summarise" %in% commandArgs(trailingOnly = TRUE)) {
  for (script in summarise_scripts) run_script(script)
}
for (script in plot_scripts) run_script(script)
