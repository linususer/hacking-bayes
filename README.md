# What is this repository about?



# Directory Structure:
```
hacking-bayes
.
├── figures
└── scripts
    ├── calculations
    └── plots
```

The simulation data of the main results of the Optional Stopping simulations with Normal and Cauchy Prior can be found under [Releases](https://github.com/linususer/hacking-bayes/releases/download/sim_data1.0/hacking-bayes-factors.duckdb).

The database is structured in four tables:

```
bf_decision_threshold
application_example
cauchy_prior
cauchy_sym
```




# Dependencies

`duckdb`
`data.table`
`BayesFactor`
