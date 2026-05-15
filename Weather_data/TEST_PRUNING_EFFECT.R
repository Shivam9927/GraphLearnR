# ============================================================================
# Test: Pruning Effect on Metrics
# ============================================================================

library(GraphLearnR)

data_path <- getwd()
w1 <- 0.0001
w2 <- 0.99

# ============================================================================
# LOAD DATA
# ============================================================================

iteration2 <- read.csv(file.path(data_path, "NOAA_iteration2.csv"), row.names=1)
dist_matrix <- read.csv(file.path(data_path, "NOAA_dist_increased.csv"), row.names=1)

N <- nrow(iteration2)

# ============================================================================
# PREPARE DATA
# ============================================================================

tmax_cols <- grep("^TMAX_", names(iteration2), value=TRUE)
TMAX <- as.matrix(iteration2[, tmax_cols])
n <- ncol(TMAX)

elevation <- iteration2$ELEVATION
tmin_cols <- grep("^TMIN_", names(iteration2), value=TRUE)
TMIN <- as.matrix(iteration2[, tmin_cols])

P1 <- matrix(rep(elevation, n), nrow=N, ncol=n)
P2 <- TMIN
P3 <- matrix(1, nrow=N, ncol=n)
P <- list(P1, P2, P3)

# ============================================================================
# CREATE GROUND TRUTH
# ============================================================================

dist_mat <- as.matrix(dist_matrix)
L_prime <- matrix(999, N, N)
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
  
  list(precision=precision, recall=recall, f_score=f_score, TP=TP, FP=FP, FN=FN)
}

# ============================================================================
# RUN GLREG
# ============================================================================

cat("\nRunning GLReg...\n")

result <- gsp_regression(
  X_noisy = TMAX,
  P = P,
  w1 = w1,
  w2 = w2,
  max_iter = 100,
  threshold = 1e-3
)

L_learned <- result$L

# ============================================================================
# TEST DIFFERENT PRUNING THRESHOLDS
# ============================================================================

cat("\n")
cat("====================================================================\n")
cat(" Testing Different Pruning Thresholds\n")
cat("====================================================================\n")

thresholds <- c(0, 0.0005, 0.001, 0.005, 0.01, 0.05)

for (thresh in thresholds) {
  L_test <- L_learned
  L_test[abs(L_test) < thresh] <- 0
  
  n_edges <- sum(abs(L_test[lower.tri(L_test)]) > 0)
  metrics <- metricsprf_partial(L_prime, L_test)
  
  cat(sprintf("\nThreshold: %.4f\n", thresh))
  cat(sprintf("  Edges: %d\n", n_edges))
  cat(sprintf("  Precision: %.4f\n", metrics$precision))
  cat(sprintf("  Recall:    %.4f\n", metrics$recall))
  cat(sprintf("  F1-Score:  %.4f\n", metrics$f_score))
}

cat("\n")
cat("====================================================================\n")
cat(" Christopher's Approach\n")
cat("====================================================================\n")
cat("\nHe evaluates on L pruned with threshold=0.0005\n")
cat("Check which threshold gives the best F1-Score above.\n\n")
