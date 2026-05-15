# ============================================================================
# NOAA Complete Analysis - All Iterations
# ============================================================================

library(GraphLearnR)

# ============================================================================
# CONFIGURATION
# ============================================================================

data_path <- getwd()
w1 <- 0.0001
w2 <- 0.99

# ============================================================================
# LOAD DISTANCE MATRIX
# ============================================================================

dist_matrix <- read.csv(
  file.path("C:/Users/sgpawar/Documents/GraphLearn R/GL_weather-main",
            "NOAA_dist_increased.csv"),
  row.names = 1
)

dist_mat <- as.matrix(dist_matrix)

# ============================================================================
# PREPARE DATA FUNCTION
# ============================================================================

prepare_data <- function(data_file) {
  data <- read.csv(file.path(data_path, data_file), row.names=1)
  
  N <- nrow(data)
  
  # Extract TMAX (signal)
  tmax_cols <- grep("^TMAX_", names(data), value=TRUE)
  TMAX <- as.matrix(data[, tmax_cols])
  n <- ncol(TMAX)
  
  # Extract regressors
  elevation <- data$ELEVATION
  tmin_cols <- grep("^TMIN_", names(data), value=TRUE)
  TMIN <- as.matrix(data[, tmin_cols])
  
  # Create regressor list
  P1 <- matrix(rep(elevation, n), nrow=N, ncol=n)
  P2 <- TMIN
  P3 <- matrix(1, nrow=N, ncol=n)
  
  list(signal = TMAX, regressors = list(P1, P2, P3), N = N, n = n)
}

# ============================================================================
# CREATE PARTIAL GROUND TRUTH
# ============================================================================

create_ground_truth <- function(dist_mat) {
  N <- nrow(dist_mat)
  
  L_prime <- matrix(999, N, N)
  L_prime[dist_mat < 40000 & dist_mat > 0] <- -1
  L_prime[dist_mat > 193121] <- 0
  diag(L_prime) <- 0
  
  L_prime
}

L_prime <- create_ground_truth(dist_mat)

# ============================================================================
# EVALUATION FUNCTION
# ============================================================================

metricsprf_partial <- function(L_prime, L_learned) {
  L_prime_lower <- L_prime[lower.tri(L_prime)]
  L_learned_lower <- L_learned[lower.tri(L_learned)]
  
  known_idx <- which(L_prime_lower != 999)
  
  L_prime_known <- L_prime_lower[known_idx]
  L_learned_known <- L_learned_lower[known_idx]
  
  true_edges <- (L_prime_known != 0)
  learned_edges <- (L_learned_known != 0)
  
  TP <- sum(true_edges & learned_edges)
  FP <- sum(!true_edges & learned_edges)
  FN <- sum(true_edges & !learned_edges)
  
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0
  recall <- if (TP + FN > 0) TP / (TP + FN) else 0
  f_score <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0
  
  list(precision=precision, recall=recall, f_score=f_score)
}

# ============================================================================
# RUN ON ALL ITERATIONS
# ============================================================================

iteration_files <- c(
  "NOAA_iteration2.csv",
  "NOAA_iteration3.csv",
  "NOAA_iteration4.csv",
  "NOAA_iteration5.csv"
)

results_list <- list()

cat("\n")
cat("====================================================================\n")
cat(" Running GLReg on All Iterations\n")
cat("====================================================================\n")

for (i in seq_along(iteration_files)) {
  cat(sprintf("\nIteration %d/%d: %s\n", i, length(iteration_files), iteration_files[i]))
  cat("--------------------------------------------------------------------\n")
  
  # Prepare data
  data <- prepare_data(iteration_files[i])
  
  # Run GLReg
  result <- gsp_regression(
    X_noisy = data$signal,
    P = data$regressors,
    w1 = w1,
    w2 = w2,
    max_iter = 100,
    threshold = 1e-3
  )
  
  L_learned <- result$L
  
  # Evaluate
  metrics <- metricsprf_partial(L_prime, L_learned)
  
  # Store results
  results_list[[i]] <- list(
    L = L_learned,
    metrics = metrics
  )
  
  cat(sprintf("Precision: %.4f\n", metrics$precision))
  cat(sprintf("Recall:    %.4f\n", metrics$recall))
  cat(sprintf("F1-Score:  %.4f\n", metrics$f_score))
}

# ============================================================================
# POOL RESULTS
# ============================================================================

cat("\n")
cat("====================================================================\n")
cat(" Pooled Results (Mean ± SD)\n")
cat("====================================================================\n")

precision_vals <- sapply(results_list, function(x) x$metrics$precision)
recall_vals <- sapply(results_list, function(x) x$metrics$recall)
f1_vals <- sapply(results_list, function(x) x$metrics$f_score)

precision_mean <- mean(precision_vals)
precision_sd <- sd(precision_vals)
recall_mean <- mean(recall_vals)
recall_sd <- sd(recall_vals)
f1_mean <- mean(f1_vals)
f1_sd <- sd(f1_vals)

cat(sprintf("\nPrecision: %.4f ± %.4f\n", precision_mean, precision_sd))
cat(sprintf("Recall:    %.4f ± %.4f\n", recall_mean, recall_sd))
cat(sprintf("F1-Score:  %.4f ± %.4f\n", f1_mean, f1_sd))

# ============================================================================
# SAVE RESULTS
# ============================================================================

pooled_results <- data.frame(
  Metric = c("Precision", "Recall", "F1-Score"),
  Mean = c(precision_mean, recall_mean, f1_mean),
  SD = c(precision_sd, recall_sd, f1_sd)
)

imputation_results <- data.frame(
  Imputation = 1:4,
  Precision = precision_vals,
  Recall = recall_vals,
  F1_Score = f1_vals
)

write.csv(pooled_results, "NOAA_pooled_results.csv", row.names=FALSE)
write.csv(imputation_results, "NOAA_all_imputations.csv", row.names=FALSE)

cat("\n✓ Saved: NOAA_pooled_results.csv\n")
cat("✓ Saved: NOAA_all_imputations.csv\n\n")
