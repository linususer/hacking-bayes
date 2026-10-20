# Simulations for the e-value study. Writes the summaries to student-research/summaries/;
# the plots are made from them by student-research/plots/e-values-plots.R.

# Set working directory
setwd(".")

# Clear workspace
rm(list = ls())
gc()
library(safestats)
# set.seed(1)

# parameter configuration
alpha <- 0.05
beta <- 0.2
deltaMin <- 9 / (sqrt(2) * 15) # minimal effect size we want to be able to detect

# check if density and T-test functions are the same
designObj <- designSafeT(
  deltaMin = deltaMin,
  alpha = alpha,
  beta = beta,
  alternative = "greater",
  testType = "paired"
)

safeTTestDensity_paired <- function(x, y, designObj, na.rm = FALSE) {
  n <- length(x)
  deltaS <- designObj[["parameter"]]
  nEff <- n
  nu <- n - 1
  meanObs <- mean(x - y, na.rm = na.rm)
  sdObs <- sd(x - y, na.rm = na.rm)
  h0 <- designObj[["h0"]]
  t <- sqrt(nEff) * (meanObs - h0)/sdObs
  logEValue <- log(dt(t, df = nu, ncp = sqrt(nEff) * deltaS)) - log(dt(t, df = nu, ncp = 0))
  #print(t)
  exp(logEValue[["mu"]])
}
# print(replicate(100, {n <- 10
# x1 <- rnorm(n, mean = 0) # true effect = deltaMin
# x2 <- rnorm(n, mean = 0)
# eVal <- safeTTestDensity_paired(x1, x2, designObj)}))

#  x1 <- rnorm(2, mean = trueDelta) # veränderbar
#  x2 <- rnorm(2, mean = 0) # unverändert
#  safeTTest(x1, x2, designObj = designObj, paired = TRUE)

# eVal_t_density <- safeTTestDensity_paired(
#   t           = tStat,
#   deltaS      = deltaMin,
#   n           = n
# )

#eVal_safestats <- safeTTest(x1, x2, designObj = designObj, paired = TRUE)$eValue

#safeTTestStat(t=3, designObj[["parameter"]], n1 = 300, n2= 300, alternative = "greater", paired = TRUE)

#eVal_t_density
#eVal_safestats

# Simulate Optional Stopping
eValOptionalStopping <- function(
    designObj, paired = TRUE,
    two_thresholded = FALSE) {
  # only paired
  if (!paired) {
    stop("Optional Stopping not implemented for non paired yet")
  }

  alpha <- designObj[["alpha"]]
  n_max <- 500 # designObj[["nPlan"]][[1]]
  eCrit1 <- 1 / alpha
  eCrit0 <- if (two_thresholded) alpha else -1 # beta => type II?

  eValues <- rep(NA_real_, n_max)
  stopped <- FALSE
  stopCount <- NA_integer_
  eValueAtStop <- NA_real_
  decision <- "no decision reached by n_max"
  
  x1 <- rnorm(2, mean = 0) # H1 true: mean = deltaMin)
  x2 <- rnorm(2, mean = 0)
  count <- 2
  while (count <= n_max) {
    eVal <- safeTTest(
      x = x1, y = x2,
      designObj = designObj, paired = TRUE
    )$eValue

    eValueAtStop <- eVal
    eValues[count] <- eVal

    if (eVal >= eCrit1) {
      stopped <- TRUE
      stopCount <- count
      decision <- "H1"
      break
    }

    if (two_thresholded && eVal <= eCrit0) {
      stopped <- TRUE
      stopCount <- count
      decision <- "H0"
      break
    }
    x2 <- c(x2, rnorm(1, mean = 0))
    x1 <- c(x1, rnorm(1, mean = 9 / (sqrt(2) * 15))) # H1 true: mean = deltaMin
    count <- count + 1
  }
  tStat <- safeTTest(
    x = x1, y = x2,
    designObj = designObj, paired = TRUE
  )
  # print(tStat)

  list(
    eValues = eValues,
    stopped = stopped,
    stopCount = stopCount,
    eValueAtStop = eValueAtStop,
    decision = decision,
    nUsed = if (stopped) stopCount else n_max
  )
}
n_max <- 1000

res <- eValOptionalStopping(n_max,
  designObj = designObj,
  two_thresholded = TRUE
)

#########################################
# replicate Optional Stopping 500 times #
#########################################

eValues <- replicate(
  n = 1000,
  expr = eValOptionalStopping(designObj, two_thresholded = TRUE)$eValueAtStop
)

print(paste0("H1: ", mean(eValues >= 20)))
print(paste0("H0: ", mean(eValues <= 0.05)))
print(paste0("ind: ", 1 - mean(eValues >= 20) - mean(eValues <= 0.05)))

# Simulate effect sizess
deltaMins <- seq(0.05, 1, 0.05)  # Start from 0.05 since deltaMin=0 causes design issues
n_reps <- 1000
# alpha <- 0.9

# Initialize results storage
all_results <- list()
summary_rows <- list()

for (i in seq_along(deltaMins)) {
  deltaMin <- deltaMins[i]
  cat(sprintf("Progress: %d/%d (deltaMin=%.2f)\n", i, length(deltaMins), deltaMin))
  designObj <- tryCatch({
    designSafeT(
      deltaMin = deltaMin,
      alpha = alpha,
      beta = beta,
      alternative = "greater",
      testType = "paired"
    )
  }, error = function(e) {
    cat("Error for deltaMin=", deltaMin, ": ", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(designObj)) {
    next
  }
  
  # Run replications and store full results
  replication_results <- replicate(
    n = n_reps,
    expr = eValOptionalStopping(designObj, two_thresholded = TRUE),
    simplify = FALSE
  )
  
  # Extract eValues at stop
  eValues <- sapply(replication_results, function(x) x$eValueAtStop)
  nUsed <- sapply(replication_results, function(x) x$nUsed)
  decisions <- sapply(replication_results, function(x) x$decision)
  
  # Calculate summary statistics
  h1_rate <- mean(eValues >= 1/alpha)
  h0_rate <- mean(eValues <= alpha)
  ind_rate <- 1 - h1_rate - h0_rate
  mean_eVal <- mean(eValues, na.rm = TRUE)
  median_eVal <- median(eValues, na.rm = TRUE)
  sd_eVal <- sd(eValues, na.rm = TRUE)
  mean_nUsed <- mean(nUsed, na.rm = TRUE)
  h1_n <- sum(decisions == "H1")
  h0_n <- sum(decisions == "H0")
  
  # Store summary row as a list for tab-separated output
  summary_rows[[as.character(deltaMin)]] <- data.frame(
    deltaMin = deltaMin,
    alpha = alpha,
    beta = beta,
    n_max = n_max,
    n_reps = n_reps,
    h1_rate = h1_rate,
    h0_rate = h0_rate,
    ind_rate = ind_rate,
    mean_eValue = mean_eVal,
    median_eValue = median_eVal,
    sd_eValue = sd_eVal,
    mean_nUsed = mean_nUsed,
    h1_n = h1_n,
    h0_n = h0_n,
    stringsAsFactors = FALSE
  )
  
  # Store full results for this deltaMin
  all_results[[as.character(deltaMin)]] <- list(
    deltaMin = deltaMin,
    designObj = designObj,
    eValues = eValues,
    nUsed = nUsed,
    decisions = decisions,
    replication_results = replication_results,
    summary = list(
      h1_rate = h1_rate,
      h0_rate = h0_rate,
      ind_rate = ind_rate,
      mean_eVal = mean_eVal,
      median_eVal = median_eVal,
      sd_eVal = sd_eVal,
      mean_nUsed = mean_nUsed,
      h1_n = h1_n,
      h0_n = h0_n
    )
  )
}

# Write summary to file in tab-separated format (matching other summaries)
summary_path <- "student-research/summaries/e-values-summary.txt"
summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
cat("Summary written to:", summary_path, "\n")

# Save data as R data object
data_path <- "data/e-values-simulation.RData"
save(all_results, alpha, beta, n_max, deltaMins, file = data_path)
cat("Data saved to:", data_path, "\n")

# Print summary to console as well
print(summary_df)

# Ensure summaries directory exists relative to hacking-bayes
dir.create("student-research/summaries", showWarnings = FALSE, recursive = TRUE)

##################################################
######## Simulation for fixed delta test design ##
##################################################


# parameter configuration
alpha <- 0.05
beta <- 0.2
fixed_deltaMin <- 9 / (sqrt(2) * 15) # minimal effect size we want to be able to detect

designObj_fixed <- designSafeT(
  deltaMin = fixed_deltaMin,
  alpha = alpha,
  beta = beta,
  alternative = "greater",
  testType = "paired"  
)

# Ensure summaries directory exists
dir.create("student-research/summaries", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Extended Optional Stopping (Incremental Data Generation)
# ==============================================================================
eValOptionalStopping <- function(
    designObj, 
    trueDelta = 0,         
    paired = TRUE, 
    two_thresholded = TRUE,
    n_max = 500,
    eCrit_upper = NULL,
    eCrit_lower = NULL) {
  
  if (!paired) stop("Optional Stopping not implemented for non-paired yet")

  # Define Thresholds
  if(is.null(eCrit_upper)) eCrit_upper <- 1 / designObj[["alpha"]]
  if(is.null(eCrit_lower)) eCrit_lower <- if(two_thresholded) designObj[["alpha"]] else -1
  
  eValues <- rep(NA_real_, n_max)
  stopped <- FALSE
  stopCount <- n_max
  opt_stop_decision <- "indecisive"
  
  # Start incrementally
  x1 <- rnorm(2, mean = trueDelta)
  x2 <- rnorm(2, mean = 0)
  count <- 2
  
  while (count <= n_max) {
    eVal <- safeTTest(
      x = x1, y = x2,
      designObj = designObj, paired = TRUE
    )$eValue
    
    eValues[count] <- eVal
    
    # Check if a threshold is crossed
    if (eVal >= eCrit_upper) {
      stopped <- TRUE
      stopCount <- count
      opt_stop_decision <- "H1"
    } else if (two_thresholded && eVal <= eCrit_lower) {
      stopped <- TRUE
      stopCount <- count
      opt_stop_decision <- "H0"
    }
    
    # Option: If a decision is reached, add the rest of the data in one go
    if (stopped) {
      remaining <- n_max - count
      if (remaining > 0) {
        x2 <- c(x2, rnorm(remaining, mean = 0))
        x1 <- c(x1, rnorm(remaining, mean = trueDelta))
        
        eValues[n_max] <- safeTTest(
          x = x1, y = x2,
          designObj = designObj, paired = TRUE
        )$eValue
      }
      break
    }
    
    # Otherwise, add one point incrementally and continue
    x2 <- c(x2, rnorm(1, mean = 0))
    x1 <- c(x1, rnorm(1, mean = trueDelta))
    count <- count + 1
  }
  
  # Final n_max Decision evaluation
  final_eVal <- eValues[n_max]
  if (final_eVal >= eCrit_upper) {
    n_max_decision <- "H1"
  } else if (two_thresholded && final_eVal <= eCrit_lower) {
    n_max_decision <- "H0"
  } else {
    n_max_decision <- "indecisive"
  }
  
  # Determine Walk Profile
  walk_profile <- paste0("OptStop ", opt_stop_decision, ", n_max ", n_max_decision)
  #print(final_eVal)
  #print(walk_profile)
  list(
    eValues = eValues,
    stopCount = stopCount,
    opt_stop_decision = opt_stop_decision,
    n_max_decision = n_max_decision,
    walk_profile = walk_profile
  )
}

# ==============================================================================
# 2. Tracked Walks Summarization (Grouped by Walk Profile)
# ==============================================================================
summarize_walks <- function(replication_results, trueDelta) {
  profiles <- sapply(replication_results, function(x) x$walk_profile)
  optDec   <- sapply(replication_results, function(x) x$opt_stop_decision)
  fixDec   <- sapply(replication_results, function(x) x$n_max_decision)
  eVal_mat <- do.call(cbind, lapply(replication_results, function(x) x$eValues))

  n_total <- length(profiles)
  n_max   <- nrow(eVal_mat)

  # the two numbers that get compared, plus who rejects that the other does not
  p_optstop_h1 <- mean(optDec == "H1")
  p_fixed_h1   <- mean(fixDec == "H1")
  p_opt_only   <- mean(optDec == "H1" & fixDec != "H1")
  p_fix_only   <- mean(optDec != "H1" & fixDec == "H1")

  summary_rows <- list()

  for (prof in unique(profiles)) {
    prof_cols <- eVal_mat[, profiles == prof, drop = FALSE]
    prof_pct  <- ncol(prof_cols) / n_total

    for (n_idx in 2:nrow(prof_cols)) {
      vals <- prof_cols[n_idx, ]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) next

      qs <- quantile(log(vals), probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        trueDelta = trueDelta,
        n_max = n_max,
        walk_profile = prof,
        profile_pct = prof_pct,
        n = n_idx,
        count_at_n = length(vals),
        log_q25 = qs[1],
        log_median = qs[2],
        log_q75 = qs[3],
        p_optstop_h1 = p_optstop_h1,
        p_fixed_h1 = p_fixed_h1,
        p_opt_only = p_opt_only,
        p_fix_only = p_fix_only,
        n_reps = n_total,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_df
}

# ==============================================================================
# Simulation 1: Extended Optional Stopping (Variable trueDelta)
# ==============================================================================
trueDeltas <- seq(0,1,0.05)
n_reps <- 1000
n_max <- 10

eCrit1 <- 1 / alpha
eCrit0 <- alpha

all_summaries <- list()

for (trueDelta in trueDeltas) {
  cat(sprintf("\n[Standard] Simulating True Delta = %.2f...\n", trueDelta))
  
  replications <- replicate(
    n = n_reps,
    expr = eValOptionalStopping(
      designObj = designObj_fixed, 
      trueDelta = trueDelta, 
      n_max = n_max,
      eCrit_upper = eCrit1,
      eCrit_lower = eCrit0
    ),
    simplify = FALSE
  )
  
  sum_df <- summarize_walks(replications, trueDelta = trueDelta)
  all_summaries[[as.character(trueDelta)]] <- sum_df
  
}

# Export combined summary
final_standard_summary <- do.call(rbind, all_summaries)
write.table(final_standard_summary, "student-research/summaries/extended-opt-stop-summary.txt", sep = "\t", row.names = FALSE)
cat("\n=> Saved Extended Optional Stopping summary to: student-research/summaries/extended-opt-stop-summary.txt\n")


# ==============================================================================
# Simulation 2: Asymmetrical Thresholds
# ==============================================================================
asym_upper <- 10        # upper bound, used by BOTH arms
asym_lower <- 0.95    # H0 bound, optional stopping only
n_maxs      <- 5:15
asym_n_reps <- 20000

all_asym_summaries <- list()

for (trueDelta in c(0)) {
  for (n_max in n_maxs){
  cat(sprintf("\n[Asymmetrical] Simulating n = %.2f...\n", n_max))

  replications_asym <- replicate(
    n = asym_n_reps,
    expr = eValOptionalStopping(
      designObj = designObj_fixed,
      trueDelta = trueDelta,
      n_max = n_max,
      eCrit_upper = asym_upper,
      eCrit_lower = asym_lower
    ),
    simplify = FALSE
  )

  sum_df_asym <- summarize_walks(replications_asym, trueDelta = trueDelta)
  all_asym_summaries[[paste(trueDelta, n_max)]] <- sum_df_asym

}
}
final_asym_summary <- do.call(rbind, all_asym_summaries)
write.table(final_asym_summary, "student-research/summaries/asym-opt-stop-summary.txt",
            sep = "\t", row.names = FALSE)
cat("\n=> Saved Asymmetrical Threshold summary to: student-research/summaries/asym-opt-stop-summary.txt\n")
# ==============================================================================
# 4. Fixed-n Simulation: decision counts for drawn sample sizes n = 2 ... 100
# ==============================================================================
# For every sample size n, draw n_reps fresh paired samples of size n,
# run safestats::safeTTest once per sample and classify the e-value:
#   H1         if eValue >= eCrit_upper (default 1 / alpha)
#   H0         if eValue <= eCrit_lower (default alpha)
#   indecisive otherwise
# The summary contains the decision counts and decision probabilities per n.
safeTTestKummer_paired <- function(x, y, designObj) {
  d <- x - y
  n <- length(d)
  nu <- n - 1
  deltaS <- unname(designObj[["parameter"]])
  t <- unname(sqrt(n) * (mean(d) - designObj[["h0"]]) / sd(d))

  a <- t^2 / (nu + t^2)
  logExpTerm <- (a - 1) * n * deltaS^2 / 2
  zArg <- -a * n * deltaS^2 / 2

  aK <- Re(hypergeo::genhypergeo(U = -nu/2, L = 1/2, zArg))
  bK <- exp(lgamma(nu/2 + 1) - lgamma((nu + 1)/2)) * sqrt(2 * n) * deltaS *
        t / sqrt(t^2 + nu) *
        Re(hypergeo::genhypergeo(U = (1 - nu)/2, L = 3/2, zArg))
  kummerSum <- aK + bK

  if (is.finite(kummerSum) && kummerSum > 0) {
    exp(logExpTerm + log(kummerSum))
  } else {
    # fallback: ratio of t densities on the log scale ("greater")
    exp(dt(t, df = nu, ncp = sqrt(n) * deltaS, log = TRUE) -
        dt(t, df = nu, log = TRUE))
  }
}


simulate_fixed_n <- function(
    designObj,
    trueDelta = 0,
    sample_sizes = 10:1000,
    n_reps = 500,
    eCrit_upper = NULL,
    eCrit_lower = NULL,
    summary_path = NULL,
    verbose = TRUE) {
  if (is.null(eCrit_upper)) eCrit_upper <- 1 / designObj[["alpha"]]
  if (is.null(eCrit_lower)) eCrit_lower <- designObj[["alpha"]]

  summary_rows <- vector("list", length(sample_sizes))

  for (i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    if (verbose) {
      cat(sprintf("[Fixed n] n=%d (%d/%d)\n", n, i, length(sample_sizes)))
    }

    decisions <- vapply(seq_len(n_reps), function(r) {
      x1 <- rnorm(n, mean = trueDelta)
      x2 <- rnorm(n, mean = 0)
      # custom density function
      
      eVal <- safeTTestKummer_paired(x1,x2,designObj = designObj)
      #eVal <- safeTTestDensity_paired(x1,x2, designObj)
      # print(eVal)
      if (eVal >= eCrit_upper) {
        "H1"
        #print("H1")
      } else if (eVal <= eCrit_lower) {
        "H0"
        #print("H0")
      } else {
        "indecisive"
        #print("Indecisive")
      }
    }, character(1))

    h1_n <- sum(decisions == "H1")
    h0_n <- sum(decisions == "H0")
    ind_n <- sum(decisions == "indecisive")

    summary_rows[[i]] <- data.frame(
      trueDelta = trueDelta,
      eCrit_upper = eCrit_upper,
      eCrit_lower = eCrit_lower,
      n = n,
      n_reps = n_reps,
      h1_n = h1_n,
      h0_n = h0_n,
      ind_n = ind_n,
      h1_prob = h1_n / n_reps,
      h0_prob = h0_n / n_reps,
      ind_prob = ind_n / n_reps,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  if (!is.null(summary_path)) {
    dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
    write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
    cat("=> Saved fixed-n summary to:", summary_path, "\n")
  }

  summary_df
}

# ==============================================================================
# Simulation 3: Fixed sample sizes (n = 2 ... 1000), 500 t-tests per n
# ==============================================================================
# --- True Delta = 0 (H0 true) ---
# set.seed(1)
cat("\n[Fixed n] Simulating True Delta = 0.00...\n")
fixed_n_summary_0 <- simulate_fixed_n(
  designObj = designObj_fixed,
  trueDelta = 0,
  sample_sizes = 10:1000,
  n_reps = 500,
  summary_path = "student-research/summaries/fixed-n-summary-trueDelta-0.00-own.txt"
)
 
# --- True Delta = 0.1 ---
cat("\n[Fixed n] Simulating True Delta = 0.10...\n")
fixed_n_summary_01 <- simulate_fixed_n(
  designObj = designObj_fixed,
  trueDelta = 0.1,
  sample_sizes = 2:1000,
  n_reps = 500,
  summary_path = "student-research/summaries/fixed-n-summary-trueDelta-0.10-own.txt"
)
 
# --- True Delta = 0.2 ---
cat("\n[Fixed n] Simulating True Delta = 0.20...\n")
fixed_n_summary_02 <- simulate_fixed_n(
  designObj = designObj_fixed,
  trueDelta = 0.2,
  sample_sizes = 2:1000,
  n_reps = 500,
  summary_path = "student-research/summaries/fixed-n-summary-trueDelta-0.20-own.txt"
)
 
# --- True Delta = 0.3 ---
cat("\n[Fixed n] Simulating True Delta = 0.30...\n")
fixed_n_summary_03 <- simulate_fixed_n(
  designObj = designObj_fixed,
  trueDelta = 0.3,
  sample_sizes = 2:1000,
  n_reps = 500,
  summary_path = "student-research/summaries/fixed-n-summary-trueDelta-0.30-own.txt"
)
 
# --- True Delta = 0.42 ---
cat("\n[Fixed n] Simulating True Delta = 0.42...\n")
fixed_n_summary_04 <- simulate_fixed_n(
  designObj = designObj_fixed,
  trueDelta = 0.42,
  sample_sizes = 10:1000,
  n_reps = 500,
  summary_path = "student-research/summaries/fixed-n-summary-trueDelta-0.42-own.txt"
)
 
# --- True Delta = 0.5 ---
cat("\n[Fixed n] Simulating True Delta = 0.50...\n")
fixed_n_summary_05 <- simulate_fixed_n(
  designObj = designObj_fixed,
  trueDelta = 0.5,
  sample_sizes = 2:1000,
  n_reps = 500,
  summary_path = "student-research/summaries/fixed-n-summary-trueDelta-0.50-own.txt"
)
 
# Export combined summary (all true deltas)
final_fixed_n_summary <- rbind(
  fixed_n_summary_0,
  fixed_n_summary_01,
  fixed_n_summary_02,
  fixed_n_summary_03,
  fixed_n_summary_04,
  fixed_n_summary_05
)
write.table(final_fixed_n_summary, "student-research/summaries/fixed-n-summary-own.txt",
            sep = "\t", quote = TRUE, row.names = FALSE)
cat("\n=> Saved combined fixed-n summary to: student-research/summaries/fixed-n-summary.txt\n")
