# clear workspace
rm(list = ls())
setwd("/home/linus/git/hacking-bayes")
# Every function call corresponds to one figure in the paper

#######################
### Catch-Up-Effect ###
#######################

source("scripts/plots/catch-up-effect-plot-functions.R")
# Figure 1 - Frequentistic Optional Stopping
freq_optional_stopping()
print("Figure 1 done")
# Figure 2 - Random Walk with Boundaries for the Symmetrical Stopping Rule
bayes_optional_stopping()
print("Figure 2 done")
# Figure 3 - Idealized Setting - Point prior vs Normal prior
point_and_normal_prior(true_mean = 0, sd_1 = 1)
print("Figure 3 done")
# Figure 4 - Idealized Setting - BF01 depending on n and mu - logscale
ns <- 10^seq(0, 4, 0.001)
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
bf01_mus_plot(ns, mus, true_var = 1, prior_var = 1, logscale = TRUE, file_name = "figures/bf01-compare-mus-logscale.pdf")
print("Figure 4 done")
# Figure 5 - Realistic Setting - Point prior vs Cauchy prior
point_and_cauchy_prior()
print("Figure 5 done")
# Figure 7 - Catch-Up Effect - Minimum everywhere except \Bar{y} = 0
ns <- 1:10000
bf01_mus_plot(ns, mus = c(0), true_var = 1, prior_var = 1, file_name = "figures/bf01-minimum-everywhere-except-mu-0.pdf", ylim = c(0, 100))
print("Figure 7 done")
# Figure 8 - Catch-Up Effect - Comparing different true variances
ns <- 1:2000
vars <- c(0.25, 1, 4, 25)
bf01_var_plot(ns, mu = 0.1, vars = vars, prior_var = 1)
print("Figure 8 done")
# Figure 9 - Catch-Up Effect - Comparing different alternative prior variances
ns <- 1:500
prior_vars <- c(0.01, 0.25, 1, 4, 25)
bf01_prior_var_plot(ns, mu = 0.1, true_var = 1, prior_vars = prior_vars)
print("Figure 9 done")
# Figure 10 - Catch-Up-Effect - Comparing different true means
ns <- 1:100
mus <- c(0.1, 0.2, 0.5, 0.8, 1)
bf01_mus_plot(ns, mus, true_var = 1, prior_var = 1)
print("Figure 10 done")
# Figure 11 - Catch-Up-Effect - normal prior vs point prior and BF10 (\sigma_1 < \sqrt{\Bar{y}^2 - \sigma^2})
point_and_normal_prior(observed_mean = 2, sd_1 = 0.01, file_name = "figures/normal-prior-sd_1<sqrt(y^2-sd^2).pdf") # \sigma = 1 => \sigma_1 = 0.01 < \sqrt{2^2 - 1} = \sqrt{3}
ns <- 1:100
bf01_mus_plot(ns, mus = c(2), true_var = 1, prior_var = 0.01, file_name = "figures/bf10-sd_1<sqrt(y^2-sd^2).pdf", ylim = c(0, 3))
print("Figure 11 done")
# Figure 12 - Catch-Up-Effect - normal prior vs point prior and BF10 (\sigma_1 > \Bar{y})
point_and_normal_prior(observed_mean = 0.5, sd_1 = 2, file_name = "figures/normal-prior-sd_1>mean.pdf") # \sigma = 1 => \sigma_1 = 1 > 0.5
ns <- 1:100
bf01_mus_plot(ns, mus = c(0.5), true_var = 1, prior_var = 1, file_name = "figures/bf10-sd_1>mean.pdf", ylim = c(0, 3))
print("Figure 12 done")
# Figure 13 - Catch-Up-Effect - Intersection points for fixed Bayes Factor thresholds BF_crit (for H0) = BF01
ns <- 1:500
bf01_intersection_points_bf_crit(ns, mu = 0.1, true_var = 1, prior_var = 1, bf_crits = c(3, 6, 10))
bf01_intersection_points_mu(ns, mus = c(0.1, 0.2, 0.5, 0.8, 1), true_var = 1, prior_var = 1, bf_crit = 3)
print("Figure 13 done")
#########################
### Optional Stopping ###
#########################
source("scripts/plots/optional-stopping-plot-functions.R")

# Figure 14 - Optional Stopping (Tendeiros Setting) - Subplot with 3 plots with P(H_0 | \mu) for different true means with n
#                                                     and P(H_0 | \mu) with true mean on x-axis (Asymmetrical Case)
mus <- c(0.1, 0.5, 0.8)
asym <- readRDS("data/from_1_asym_simulation_data")
sim_histograms(asym, mus, xlim_max = 60, ylimh0 = 2500, ylimh1 = 2500)
decision_prob_curve(case = "asymmetrical", n_start = 1, approx = FALSE)
print("Figure 14 done")
# Figure 15 - Optional Stopping (Tendeiros Setting) - Subplot with 3 plots with P(H_0 | \mu) for different true means with n
#                                                     and P(H_0 | \mu) with true mean on x-axis (Symmetrical Case).
mus <- c(0.1, 0.5, 0.8)
sym <- readRDS("data/from_1_sym_simulation_data")
sim_histograms(sym, mus, xlim_max = 60, ylimh0 = 2500, ylimh1 = 2500)
decision_prob_curve(case = "symmetrical", n_start = 1, approx = FALSE)
print("Figure 15 done")
# Figure 16 - Optional Stopping (Tendeiros Setting) - P(H_0 | \mu) with true mean on x-axis with both stopping rules
decision_prob_curve(case = "both", n_start = 1, approx = FALSE)
print("Figure 16 done")
# Figure 17 - Optional Stopping Approximation - Subplot: Critical Means for H0/H1 and 95% CI around \mu = 0.1 with normal distributed means around \mu = 0.1
#                                               and random walk with critical means for decision making
ns <- 1:150
critical_means_boundary(ns, mu = 0.4, prior_var = 1, bf_crit = 3, step = 10, symmetric = TRUE)
critical_means_boundary(ns, mu = 0.4, prior_var = 1, bf_crit = 3, step = 10, symmetric = FALSE)
print("Figure 17 done")
# Figure 18 - Optional Stopping Approximation - Decision probability for asymmetrical rule
decision_prob_curve(case = "asymmetrical", n_start = 1, approx = TRUE)
print("Figure 18 done")
# Figure 19 - Optional Stopping Approximation - Decision probability for symmetrical rule
decision_prob_curve(case = "symmetrical", n_start = 1, approx = TRUE)
print("Figure 19 done")
# Figure 20 - Optional Stopping Approximation - Subplot: Optimize symmetrical rule for different k
#                                               and plot decision probability for H0 with chosen k
print("This will take some time ~ 30 minutes on my machine (average laptop)")
sim_diff_n_start()
decision_prob_curve(case = "symmetrical", n_start = 6, approx = TRUE)
print("Figure 20 done")
######################
# Realistic Setting #
######################
source("scripts/plots/realistic-plot-functions.R")
# Figure 21 - Optional Stopping (Realistic Setting) - Subplot with 3 plots with P(H_0 | \mu) for different true means with n
#                                                     and P(H_0 | \mu) with true mean on x-axis (Symmetrical Case).
mus <- c(0.1, 0.5, 0.8)
realistic_sim_histograms(bf_crit = 3, mus, r_val = 1 / sqrt(2))
realistic_sim_decision_prob_curve(bf_crit = 3, r_val = 1 / sqrt(2))
print("Figure 21 done")
# Figure 22 - Optional Stopping (Realistic Setting) - Show for different r and BF_crit the decision probability for H0 given mu
realistic_sim_overview_plot()
print("Figure 22 done")

source("scripts/plots/catch-up-effect-plot-functions.R")
# Figure 23 - Application Example
ns <- 1:1000
bf01_mus_plot(ns, mus = c(0.1), true_var = 2^2, prior_var = 3^2, file_name = "figures/bf01-application-h0-max.pdf", ylim = c(0, 30))
bf01_mus_plot(ns, mus = c(0.1), true_var = 1^2, prior_var = 0.15^2, file_name = "figures/bf01-application-h1-max.pdf", ylim = c(0, 2))
source("scripts/plots/optional-stopping-plot-functions.R")
h1max <- readRDS("data/application_sim_data_from1_to_Inf")
sim_histograms(df = h1max, mu = 0.1, xlim_max = 60, ylimh0 = 50, ylimh1 = 5000)
h0max <- readRDS("data/application_sim_data_from4_to_2590")
sim_histograms(h0max, 0.1, xlim_max = 60, ylimh0 = 6000, ylimh1 = 2000)
print("Figure 23 done")

# Figure 6 - Realistic Setting - Simulation for Cauchy prior to get BFs
# #                                depending on n and mu
# cauchy_plot()
# print("Figure 6 done")
