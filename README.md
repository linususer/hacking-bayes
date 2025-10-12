# What is this repository about?

Here are simulation data and plot functions to compare Optional Stopping and fixed sample size designs. Both, a bachelor thesis and a research project are build upon this. 

# Directory and DuckDB structure:

`T` denotes files that are exclusive to the Bachelor thesis, `P` denotes files that are exclusively done for a research project. `P` builds up on results of `T`, therefore these are bundled in this directory together. 

```
hacking-bayes
.
├── figures
└── scripts
    ├── calculations
        ├── applications-sim.R [T]
            -> generates application examples
        ├── cauchy-function-sim.R [T]
            -> generates random walks for cauchy prior
        ├── realisic-fix-sim.R [P] 
            -> generates the main Bayesian t-test data for fixed sample size tests
        ├── realistic-opt-sim-par.R [T]
            -> generates the main Bayesian t-test data for optional stopping tests
        ├── rouder-simulation.R [T]
            -> generates the main idealised setting data for optional stopping tests
        ├── sanborn_prob_replication.R [P]
            -> generates replication of Sanborn et al. (2014)
        └── sanborn_replication.R [P]
            -> generates replication of Sanborn et al. (2014)
    └── plots
        ├── binom-example.R [T]
            -> sanborn example plot for binomial distribution
        ├── catch-up-effect-plot-functions.R [T]
            -> Catch Up Effect plots
        ├── cauchy-plots.R [T]
            -> Plots for Cauchy prior related stuff
        ├── define_colors.R [T,P]
            -> Color definitions
        ├── fixed-plots.R [P]
            -> fixed simulation results visualised in comparison to Optional Stopping results
        ├── optional-stopping-plot-functions.R [T]
            -> different optional stopping visualisation functions
        ├── plot-all-figures.R [T]
            -> DEPRECATED
        ├── presentation.R [T]
            -> some visualisations for a presentation
        ├── random_walk_priors.R
            -> Visualisation for random walk and point priors
        ├── realistic-plot-functions.R [T]
            -> visualisations for Optional Stopping used in a poster at TEAP 2025 (see https://osf.io/yx8ng/files/f9bm3)
        ├── sanborn-replication-plots.R [P]
            -> replication plots for Sanborn et al. (2014)
        └── steele_replication.py [P]
            -> Visualisation tree for the specific coinflip example for Steele (2013)
```

The simulation data of the main results of the simulations is provided under <linus-szillat.de/ressources/hacking-bayes.duckdb>.

The database is structured in these tables:

```
├── application_example [T]
├── bf_decision_threshold [T]
├── cauchy_prior [T]
├── cauchy_sym [T]
├── cauchy_sym_fixed_size [P]
├── cauchy_sym_fixed_size_bf_crit [P]
├── cauchy_sym_fixed_size_r [P]
├── cauchy_sym_fixed_size_extra [P]
├── sanborn_probs_replication [P]
└── sanborn_replication [P]
```




# Main dependencies

`duckdb`
`data.table`
`BayesFactor`

For plotly visualisations Python and pip dependencies might be required.

# References
- Sanborn, A.N., Hills, T.T. The frequentist implications of optional stopping on Bayesian hypothesis tests. Psychon Bull Rev 21, 283–300 (2014). https://doi.org/10.3758/s13423-013-0518-9
- Steele, K. (2013). Persistent Experimenters, Stopping Rules, and Statistical Inference. Erkenntnis, 78 (4), 937–961. https://doi.org/10.1007/s10670-012-9388-1