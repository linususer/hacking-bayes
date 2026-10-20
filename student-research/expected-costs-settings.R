# Settings shared by student-research/summarise/expected-costs.R and plots/expected-costs.R:
# the effect size scale and the cost matrices of the scenarios.
r_vals <- c(0.5, 1, 2) / sqrt(2)

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


###########################################
############ WEIGHT COSTS #################
###########################################


# NEW: map each BF_crit to its purpose-built cost matrix
bf_crit_cost_matrices <- list(
    `3`  = cost_matrix4,
    `6`  = cost_matrix5,
    `10` = cost_matrix6
)
