# What is this repository about?

Here are simulation data and plot functions to compare Optional Stopping and fixed sample size designs. Both, a bachelor thesis and a research project are build upon this. 

# Directory and DuckDB structure:

Every analysis that uses the simulation data is split into two stages:

- `<project>/summarise/<topic>.R` queries `data/*.duckdb` with `*_summary()` functions and writes
  tab-separated summaries to `<project>/summaries/`.
- `<project>/plots/<topic>.R` draws the figures into `<project>/figures/` with `plot_*_summary()`
  functions and only reads `<project>/summaries/`, so all figures can be rebuilt without the database.

`Rscript run-all.R` rebuilds all figures from the committed summaries; `Rscript run-all.R summarise`
rebuilds the summaries from the database first. The simulations in `<project>/calculations/` are never
run by it.

The repository bundles three projects that build on each other. Every script is run from the repository root (paths such as `data/` and `shared/` are relative to it).

- **Thesis** – Bayesian Optional Stopping with the Bayesian t-test (`BayesFactor`).
- **Lab project** – fixed sample size designs compared with Optional Stopping, builds on the thesis data.
- **Student research** – e-values, expected costs and large effect sizes, builds on thesis and lab project data.

```
hacking-bayes
├── run-all.R                              -> rebuilds all figures (and with `summarise` all summaries)
├── data/                                  (not tracked, see below)
├── shared/
│   ├── define_colors.R                    -> color definitions (Okabe & Ito)
│   ├── queries.R                          -> DuckDB queries shared by lab project and student research
│   └── summary-io.R                       -> write_summary() / read_summary()
├── thesis/
│   ├── createDB.R                         -> migration of early simulation files into DuckDB
│   ├── calculations/
│   │   ├── cauchy-function-sim.R          -> random walks for the Cauchy prior
│   │   └── realistic-opt-sim-par.R        -> main Bayesian t-test data for Optional Stopping
│   ├── summarise/ and plots/
│   │   ├── realistic-sim.R                -> Optional Stopping with the Cauchy prior, poster at TeaP 2025 (see https://osf.io/yx8ng/files/f9bm3)
│   │   ├── decision-prob-bf-crit.R        -> decision probability for H0 by BF_crit and r
│   │   └── optional-stopping.R            -> histograms and decision probabilities for `bf_decision_threshold`
│   ├── plots/ (no database needed)
│   │   ├── binom-example.R                -> Sanborn example plot for the binomial distribution
│   │   └── catch-up-effect-plot-functions.R -> Catch Up Effect plots
│   ├── summaries/
│   └── figures/
├── lab-project/
│   ├── calculations/
│   │   └── realistic-fix-sim-par.R        -> main Bayesian t-test data for fixed sample size tests
│   ├── summarise/ and plots/
│   │   └── fixed-size.R                   -> fixed sample size results compared with Optional Stopping
│   ├── plots/ (no database needed)
│   │   ├── random_walk_priors.R           -> random walk and point prior visualisations
│   │   └── steele_replication.py          -> decision tree for the coin flip example of Steele (2013)
│   ├── summaries/
│   └── figures/                           (report/ holds the figures used in the report)
├── student-research/
│   ├── expected-costs-settings.R          -> cost matrices used by both stages of the expected costs
│   ├── calculations/
│   │   ├── e-values-simulation.R          -> e-value simulations (safestats), writes the summaries
│   │   └── extension_sim.R                -> fixed N and Optional Stopping for large effect sizes
│   ├── summarise/ and plots/
│   │   └── expected-costs.R               -> expected costs of fixed N vs. Optional Stopping
│   ├── plots/
│   │   └── e-values.R                     -> e-value plots from the simulation summaries
│   ├── summaries/
│   └── figures/
└── old/                                   -> results that are no longer needed
    ├── thesis/                            -> everything depending on the normal / point prior Bayes factor
    │                                         (bayes-factor-functions.R) and its figures
    └── lab-project/                       -> replication of Sanborn et al. (2014)
```

The simulation data of the main results of the simulations is provided under <https://linus-szillat.de/ressources/hacking-bayes.duckdb>.

The database is structured in these tables (`T` thesis, `P` lab project, `S` student research, `old` no longer used):

```
data/hacking-bayes.duckdb
├── application_example [old]           (created by old/thesis/application-sim.R)
├── bf_decision_threshold [T]             (created by old/thesis/rouder-simulation.R)
├── cauchy_prior [T]
├── cauchy_sym [T]                         (also used by P and S)
├── cauchy_sym_fixed_size [P]              (also used by S)
├── cauchy_sym_fixed_size_bf_crit [P]
├── cauchy_sym_fixed_size_r [P]
├── cauchy_sym_fixed_size_extra [P]
├── sanborn_probs_replication [old]
└── sanborn_replication [old]

data/large_effects.duckdb [S]
├── cauchy_sym
└── cauchy_sym_fixed_size

data/e-values-simulation.RData [S]
```

# Main dependencies

`duckdb`
`data.table`
`BayesFactor`
`safestats` (student research)

For plotly visualisations Python and pip dependencies might be required.

# References
- Sanborn, A.N., Hills, T.T. The frequentist implications of optional stopping on Bayesian hypothesis tests. Psychon Bull Rev 21, 283–300 (2014). https://doi.org/10.3758/s13423-013-0518-9
- Steele, K. (2013). Persistent Experimenters, Stopping Rules, and Statistical Inference. Erkenntnis, 78 (4), 937–961. https://doi.org/10.1007/s10670-012-9388-1
