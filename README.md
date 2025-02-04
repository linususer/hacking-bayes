# What is this repository about?



# Directory Structure:
```
hacking-bayes
.
├── data
├── figures
└── scripts
    ├── calculations
    └── plots
```

Under `data` the main results of the Optional Stopping simulations with Normal and Cauchy Prior can be found in `hacking-bayes-results.db`. 

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