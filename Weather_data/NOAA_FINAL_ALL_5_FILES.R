# ============================================================================
# NOAA Complete Analysis - All 5 MICE Imputations
# ============================================================================

library(GraphLearnR)

data_path <- getwd()

# ============================================================================
# CONFIGURATION
# ============================================================================

w1 <- 0.001
w2 <- 0.99

cat("\nParameters: w1 = %.4f, w2 = %.2f\n", w1, w2)

# ============================================================================
# LOAD DISTANCE MATRIX
# ============================================================================

dist_matrix <- read.csv(file.path(data_path, "NOAA_dist_increased.csv"), row.names=1)
dist_mat <- as.matrix(dist_matrix)

# ============================================================================
# CREATE PARTIAL GROUND TRUTH
# ============================================================================

L_prime <- matrix(999, 20, 20)
L_prime[dist_mat < 40000 & dist_mat > 0] <- -1
L_prime[dist_mat > 193121] <- 0
diag(L_prime) <- 0

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
# PREPARE DATA FUNCTION
# ============================================================================

prepare_and_run <- function(filename) {
  data <- read.csv(file.path(data_path, filename), row.names=1)
  
  N <- nrow(data)
  
  # Extract signal
  tmax_cols <- grep("^TMAX_", names(data), value=TRUE)
  TMAX <- as.matrix(data[, tmax_cols])
  n <- ncol(TMAX)
  
  # Extract regressors
  elevation <- data$ELEVATION
  tmin_cols <- grep("^TMIN_", names(data), value=TRUE)
  TMIN <- as.matrix(data[, tmin_cols])
  
  # Create regressor list (3 regressors)
  P1 <- matrix(rep(elevation, n), nrow=N, ncol=n)
  P2 <- TMIN
  P3 <- matrix(1, nrow=N, ncol=n)
  P <- list(P1, P2, P3)
  
  # Run GLReg
  result <- gsp_regression(
    X_noisy = TMAX,
    P = P,
    w1 = w1,
    w2 = w2,
    max_iter = 100,
    threshold = 1e-3
  )
  
  L_learned <- result$L
  
  # Evaluate
  metrics <- metricsprf_partial(L_prime, L_learned)
  
  list(L = L_learned, metrics = metrics)
}

# ============================================================================
# RUN ON ALL 5 FILES
# ============================================================================

iteration_files <- c(
  "NOAA_increased.csv",
  "NOAA_iteration2.csv",
  "NOAA_iteration3.csv",
  "NOAA_iteration4.csv",
  "NOAA_iteration5.csv"
)

results_list <- list()

cat("\n")
cat("====================================================================\n")
cat(" Running on All 5 MICE Imputations\n")
cat("====================================================================\n")

for (i in seq_along(iteration_files)) {
  cat(sprintf("\nImputation %d: %s\n", i, iteration_files[i]))
  cat("--------------------------------------------------------------------\n")
  
  result <- prepare_and_run(iteration_files[i])
  results_list[[i]] <- result
  
  cat(sprintf("Precision: %.4f\n", result$metrics$precision))
  cat(sprintf("Recall:    %.4f\n", result$metrics$recall))
  cat(sprintf("F1-Score:  %.4f\n", result$metrics$f_score))
}

# ============================================================================
# POOL RESULTS
# ============================================================================

precision_vals <- sapply(results_list, function(x) x$metrics$precision)
recall_vals <- sapply(results_list, function(x) x$metrics$recall)
f1_vals <- sapply(results_list, function(x) x$metrics$f_score)

precision_mean <- mean(precision_vals)
precision_sd <- sd(precision_vals)
recall_mean <- mean(recall_vals)
recall_sd <- sd(recall_vals)
f1_mean <- mean(f1_vals)
f1_sd <- sd(f1_vals)

cat("\n")
cat("====================================================================\n")
cat(" Pooled Results (Mean ± SD)\n")
cat("====================================================================\n")
cat(sprintf("\nPrecision: %.4f ± %.4f\n", precision_mean, precision_sd))
cat(sprintf("Recall:    %.4f ± %.4f\n", recall_mean, recall_sd))
cat(sprintf("F1-Score:  %.4f ± %.4f\n", f1_mean, f1_sd))

cat("\n")
cat("====================================================================\n")
cat(" Target (Christopher's Paper)\n")
cat("====================================================================\n")
cat("\nPrecision: 0.4689 ± 0.0418\n")
cat("Recall:    1.0000 ± 0.0000\n")
cat("F1-Score:  0.6382 ± 0.0418\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

pooled_results <- data.frame(
  Metric = c("Precision", "Recall", "F1-Score"),
  Mean = c(precision_mean, recall_mean, f1_mean),
  SD = c(precision_sd, recall_sd, f1_sd)
)

imputation_results <- data.frame(
  Imputation = 1:5,
  File = iteration_files,
  Precision = precision_vals,
  Recall = recall_vals,
  F1_Score = f1_vals
)

write.csv(pooled_results, "NOAA_final_pooled.csv", row.names=FALSE)
write.csv(imputation_results, "NOAA_final_all_imputations.csv", row.names=FALSE)

cat("\n✓ Saved: NOAA_final_pooled.csv\n")
cat("✓ Saved: NOAA_final_all_imputations.csv\n\n")
