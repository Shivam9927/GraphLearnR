# ============================================================================
# NOAA Analysis - Clean Test Script
# ============================================================================

library(GraphLearnR)

# ============================================================================
# CONFIGURATION
# ============================================================================

data_path <- getwd()
w1 <- 0.0001
w2 <- 0.99

# ============================================================================
# LOAD DATA
# ============================================================================

iteration2 <- read.csv(file.path(data_path, "NOAA_iteration2.csv"), row.names=1)
dist_matrix <- read.csv(file.path(data_path, "NOAA_dist_increased.csv"), row.names=1)

N <- nrow(iteration2)
n <- length(grep("^TMAX_", names(iteration2)))

# ============================================================================
# PREPARE SIGNAL AND REGRESSORS
# ============================================================================

# Extract TMAX (signal)
tmax_cols <- grep("^TMAX_", names(iteration2), value=TRUE)
TMAX <- as.matrix(iteration2[, tmax_cols])

# Extract regressors
elevation <- iteration2$ELEVATION
tmin_cols <- grep("^TMIN_", names(iteration2), value=TRUE)
TMIN <- as.matrix(iteration2[, tmin_cols])

# Create regressor list
P1 <- matrix(rep(elevation, n), nrow=N, ncol=n)
P2 <- TMIN
P3 <- matrix(1, nrow=N, ncol=n)

P <- list(P1, P2, P3)

# ============================================================================
# CREATE PARTIAL GROUND TRUTH
# ============================================================================

dist_mat <- as.matrix(dist_matrix)

L_prime <- matrix(999, N, N)
L_prime[dist_mat < 40000 & dist_mat > 0] <- -1
L_prime[dist_mat > 193121] <- 0
diag(L_prime) <- 0

n_known_edges <- sum(L_prime[lower.tri(L_prime)] == -1)
n_known_nonedges <- sum(L_prime[lower.tri(L_prime)] == 0 & dist_mat[lower.tri(dist_mat)] > 0)

cat(sprintf("Ground truth: %d edges, %d non-edges\n", n_known_edges, n_known_nonedges))

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
# EVALUATE
# ============================================================================

cat("\nEvaluation:\n")
cat("--------------------------------------------------------------------\n")

metrics <- metricsprf_partial(L_prime, L_learned)

cat(sprintf("TP: %d, FP: %d, FN: %d\n", metrics$TP, metrics$FP, metrics$FN))
cat(sprintf("Precision: %.4f\n", metrics$precision))
cat(sprintf("Recall: %.4f\n", metrics$recall))
cat(sprintf("F1-Score: %.4f\n", metrics$f_score))

# ============================================================================
# DIAGNOSTICS
# ============================================================================

cat("\nDiagnostics:\n")
cat("--------------------------------------------------------------------\n")
cat(sprintf("Trace(L): %.4f (should be %.0f)\n", sum(diag(L_learned)), N))
cat(sprintf("Edges learned: %d\n", sum(abs(L_learned[lower.tri(L_learned)]) > 1e-10)))

cat("\n")
