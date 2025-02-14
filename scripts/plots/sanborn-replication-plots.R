setwd("/home/linus/git/lab-project")
# clear workspace
rm(list = ls())
source("scripts/plots/define_colors.R")

# Plot the decision cumulative decision probability for preferred hypothesis H0
# dependent on number of trials n
plotSanbornSim <- function(df_lst, bf_crit, file, sym = FALSE) {
    pdf(file)
    plot(0, 0,
        xlim = c(0, 100), ylim = c(0, 1), type = "n",
        main = "Decision probability (Sanborn replication)",
        ylab = bquote("Cumulative Decision Probability"), xlab = "Number of Trials"
    )
    for (i in 1:length(df_lst)) {
        df <- df_lst[[i]]
        max_count_stop <- max(df["stop_count"])
        h0_data <- df[df["decision"] == 0 & df["bf_crit"] == bf_crit, ]
        h1_data <- df[df["decision"] == 1 & df["bf_crit"] == bf_crit, ]
        repetitions <- nrow(h0_data) + nrow(h1_data)
        stop_probs <- c()
        for (j in 1:max_count_stop) {
            h0_count <- nrow(h0_data[h0_data["stop_count"] == j, ])
            stop_prob <- h0_count / repetitions
            stop_probs <- c(stop_probs, stop_prob)
        }
        lines(1:max_count_stop, cumsum(stop_probs), col = my_colors[i])
        if (max_count_stop < 100) {
            lines(max_count_stop:100, rep(cumsum(stop_probs)[max_count_stop], 101 - max_count_stop), col = my_colors[i])
        }
        if (sym) {
            stop_probs <- c()
            for (i in 1:max_count_stop) {
                h1_count <- nrow(h1_data[h1_data["stop_count"] == i, ])
                stop_prob <- h1_count / repetitions
                stop_probs <- c(stop_probs, stop_prob)
            }
            lines(1:max_count_stop, cumsum(stop_probs), col = my_colors[2])
            # Extend the plot beyond max_count_stop
            if (max_count_stop < 100) {
                lines(max_count_stop:100, rep(cumsum(stop_probs)[max_count_stop], 101 - max_count_stop), col = my_colors[2])
            }
        }
        # abline(v = max_count_stop, col = my_colors[3])
    }
    legend("topright", legend = c("Decision for P0", "Decision for P1"), col = my_colors, lwd = 2)
    dev.off()
}

plotSanbornProbs <- function(df_list, file) {
    pdf(file)
    plot(0, 0,
        xlim = c(0.25, 0.75), ylim = c(0, 1), type = "n",
        main = "1 to 100,000 Trials vs. 20 to 200 Trials",
        ylab = bquote("Probability of Success"), xlab = "True Probability of Heads"
    )
    for (i in 1:2) {
        df <- df_list[[i]]
        lines(unlist(df[["p_true"]]), unlist(df["decision_prob0"]), col = my_colors[i])
    }
    abline(h = 1 / 10, col = my_colors[3], lty = 2)
    legend("topright", legend = c("1 to 100,000", "20 to 200"), col = my_colors, lwd = 2)
    dev.off()
}


#####################
##### LOAD DATA #####
#####################

# Open a connection to the duckdb
con <- dbConnect(duckdb(), "data/results.duckdb", read_only = TRUE)
# Load the data
asym_heads <- list()
asym_heads[[1]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                    WHERE p0 = 0.75
                                    AND p1 = 0.25
                                    AND p_true = 0.75
                                    AND symmetrical = FALSE
                                    AND alt_comp = FALSE")
asym_heads[[2]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                    WHERE p0 = 0.25
                                    AND p1 = 0.75
                                    AND p_true = 0.75
                                    AND symmetrical = FALSE
                                    AND alt_comp = FALSE")
asym_comp <- list()
asym_comp[[1]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                    WHERE p0 = 0.75
                                    AND p1 = 0.25
                                    AND p_true = 0.5
                                    AND symmetrical = FALSE
                                    AND alt_comp = FALSE")
asym_comp[[2]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                    WHERE p0 = 0.25
                                    AND p1 = 0.75
                                    AND p_true = 0.5
                                    AND symmetrical = FALSE
                                    AND alt_comp = FALSE")
asym_alt_comp <- list()
asym_alt_comp[[1]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                        WHERE p0 = 0.75
                                        AND p1 = 0.25
                                        AND symmetrical = FALSE
                                        AND alt_comp = TRUE
                                        AND alt_id = 'asym_comp1'")
asym_alt_comp[[2]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                        WHERE p0 = 0.25
                                        AND p1 = 0.75
                                        AND symmetrical = FALSE
                                        AND alt_comp = TRUE
                                        AND alt_id = 'asym_comp2'")
asym_all_comp <- list()
asym_all_comp[[1]] <- asym_comp[[1]]
asym_all_comp[[2]] <- asym_comp[[2]]
asym_all_comp[[3]] <- asym_alt_comp[[1]]
asym_all_comp[[4]] <- asym_alt_comp[[2]]

sym_heads <- list()
sym_heads[[1]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                WHERE p0 = 0.75
                                AND p1 = 0.25
                                AND p_true = 0.75
                                AND symmetrical = TRUE
                                AND alt_comp = FALSE")
sym_comp <- list()
sym_comp[[1]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                WHERE p0 = 0.75
                                AND p1 = 0.25
                                AND p_true = 0.5
                                AND symmetrical = TRUE
                                AND alt_comp = FALSE")
sym_alt_comp <- list()
sym_alt_comp[[1]] <- dbGetQuery(con, "SELECT * FROM sanborn_replication
                                    WHERE p0 = 0.75
                                    AND p1 = 0.25
                                    AND symmetrical = TRUE
                                    AND alt_comp = TRUE
                                    AND alt_id = 'sym_comp1'")

# Query data for sanborn_probs
sanborn_probs <- list()
sanborn_probs[[1]] <- dbGetQuery(con, " SELECT p_true,
                                        COUNT(CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS decision_prob0,
                                        COUNT(CASE WHEN decision = 1 THEN 1 END) * 1.0 / COUNT(*) AS decision_prob1,
                                        COUNT(CASE WHEN decision = 2 THEN 1 END) * 1.0 / COUNT(*) AS decision_prob2
                                    FROM sanborn_probs_replication
                                    WHERE trial_start = 1
                                    AND trial_end = 100000
                                    GROUP BY p_true
                                    ORDER BY p_true")
sanborn_probs[[2]] <- dbGetQuery(con, " SELECT p_true,
                                        COUNT(CASE WHEN decision = 0 THEN 1 END) * 1.0 / COUNT(*) AS decision_prob0,
                                        COUNT(CASE WHEN decision = 1 THEN 1 END) * 1.0 / COUNT(*) AS decision_prob1,
                                        COUNT(CASE WHEN decision = 2 THEN 1 END) * 1.0 / COUNT(*) AS decision_prob2
                                    FROM sanborn_probs_replication
                                    WHERE trial_start = 20
                                    AND trial_end = 200
                                    GROUP BY p_true
                                    ORDER BY p_true")


# Close the connection
dbDisconnect(con)

#####################
##### PLOT DATA #####
#####################
plotSanbornSim(asym_all_comp, bf_crit = 10, "figures/sanborn-asym-all-comp-decision-prob.pdf")
plotSanbornSim(asym_heads, bf_crit = 10, "figures/sanborn-asym-heads-decision-prob.pdf")
plotSanbornSim(asym_comp, bf_crit = 10, "figures/sanborn-asym-comp-decision-prob.pdf")
plotSanbornSim(asym_alt_comp, bf_crit = 10, "figures/sanborn-asym-alt-comp-decision-prob.pdf")
plotSanbornSim(sym_heads, bf_crit = 10, "figures/sanborn-sym-heads-decision-prob.pdf", sym = TRUE)
plotSanbornSim(sym_comp, bf_crit = 10, "figures/sanborn-sym-comp-decision-prob.pdf", sym = TRUE)
plotSanbornSim(sym_alt_comp, bf_crit = 10, "figures/sanborn-sym-alt-comp-decision-prob.pdf", sym = TRUE)

plotSanbornProbs(sanborn_probs, "figures/sanborn-probs.pdf")
