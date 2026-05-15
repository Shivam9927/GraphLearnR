# ============================================================================
# NOAA Python-matched R implementation
# Mirrors the custom Python approach in Final_NOAA_Code.ipynb / Estimating_Laplacian
# ============================================================================

suppressPackageStartupMessages({
  library(CVXR)
})

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

data_path <- getwd()
iteration_files <- c(
  "NOAA_increased.csv",
  "NOAA_iteration2.csv",
  "NOAA_iteration3.csv",
  "NOAA_iteration4.csv",
  "NOAA_iteration5.csv"
)

w1 <- 0.001
w2 <- 0.99
max_iter <- 100
prune_threshold <- 5e-4   # matches Python custom GSPRegression default
convergence_tol <- 1e-4   # matches Python stopping rule

# ----------------------------------------------------------------------------
# Helpers that mirror the Python code
# ----------------------------------------------------------------------------

lowertrian <- function(Lap) {
  Lap[lower.tri(Lap, diag = FALSE)]
}

vecF <- function(S) {
  as.vector(as.matrix(S))
}

Pbeta <- function(P, beta) {
  # Mirrors the later Python definition that stacks each regressor in Fortran order
  P_stack <- do.call(cbind, lapply(P, function(Pi) vecF(Pi)))
  beta <- as.numeric(beta)
  as.vector(P_stack %*% beta)
}

pstack <- function(P) {
  do.call(cbind, lapply(P, function(p) vecF(p)))
}

prune <- function(L, threshold) {
  temp <- L
  temp[abs(temp) < threshold] <- 0
  temp
}

normalized_mutual_info_binary <- function(x, y) {
  x <- as.integer(x)
  y <- as.integer(y)

  tab_xy <- table(factor(x, levels = c(0, 1)), factor(y, levels = c(0, 1)))
  n <- sum(tab_xy)
  if (n == 0) return(NA_real_)

  pxy <- tab_xy / n
  px <- rowSums(pxy)
  py <- colSums(pxy)

  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(0)

  mi <- 0
  for (i in seq_len(nrow(pxy))) {
    for (j in seq_len(ncol(pxy))) {
      if (pxy[i, j] > 0) {
        mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }

  mi / sqrt(hx * hy)
}

metricsprf_partial <- function(Lprime, L) {
  # Mirrors the later NOAA notebook metric version with 999 as unknown
  Lprime_tmp <- lowertrian(Lprime)
  L_tmp <- lowertrian(L)

  known_indices <- which(Lprime_tmp != 999)
  if (length(known_indices) == 0) {
    return(list(precision = -1, recall = -1, f1 = -1, nmi = -1, accuracy = -1,
                tp = NA_integer_, fp = NA_integer_, fn = NA_integer_, tn = NA_integer_))
  }

  Lprime_known <- Lprime_tmp[known_indices]
  L_known <- L_tmp[known_indices]

  true_edges <- Lprime_known != 0
  learned_edges <- L_known != 0

  if (sum(learned_edges) == 0) {
    return(list(precision = -1, recall = -1, f1 = -1, nmi = -1, accuracy = -1,
                tp = 0L, fp = 0L, fn = sum(true_edges), tn = sum(!true_edges)))
  }

  tp <- sum(true_edges & learned_edges)
  fp <- sum(!true_edges & learned_edges)
  fn <- sum(true_edges & !learned_edges)
  tn <- sum(!true_edges & !learned_edges)

  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  nmi <- normalized_mutual_info_binary(true_edges, learned_edges)

  list(
    precision = precision,
    recall = recall,
    f1 = f1,
    nmi = nmi,
    accuracy = accuracy,
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn
  )
}

getLhatprime <- function(Lhat_in, Lprime_in) {
  Lhatprime <- Lhat_in
  known_mask <- Lprime_in != 999
  Lhatprime[known_mask] <- Lprime_in[known_mask]
  Lhatprime
}

# ----------------------------------------------------------------------------
# Optimization steps translated from Python
# ----------------------------------------------------------------------------

optimize_L <- function(hatS, P, beta, w1, w2, solver = "SCS") {
  N <- nrow(hatS)
  n <- ncol(hatS)
  
  L <- Variable(c(N, N), symmetric = TRUE)
  
  R <- hatS - P[[1]] * beta[1] - P[[2]] * beta[2] - P[[3]] * beta[3]
  one_vec <- matrix(1, nrow = N, ncol = 1)
  off_diag_mask <- matrix(1, N, N) - diag(N)
  
  constraints <- list(
    sum(diag(L)) == N,
    multiply(off_diag_mask, L) <= 0,
    diag(L) >= 0,
    L %*% one_vec == 0
  )
  
  S_mat <- R %*% t(R)
  
  objective <- Minimize(
    w1 * sum_entries(L * S_mat) +
      w2 * sum_squares(L)
  )
  
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = solver, verbose = FALSE)
  
  status_now <- result$status
  if (status_now %in% c("infeasible", "unbounded")) {
    stop(sprintf("optimize_L failed with solver status: %s", status_now))
  }
  
  out <- result$getValue(L)
  out <- (out + t(out)) / 2
  out
}

optimize_hatS <- function(L, P, beta, S, w1) {
  N <- nrow(S)
  n <- ncol(S)

  InXL <- kronecker(diag(n), L)
  term1 <- w1 * InXL + diag(N * n)
  term2 <- MASS::ginv(term1)
  term3 <- vecF(S) + w1 * (InXL %*% Pbeta(P, beta))

  matrix(term2 %*% term3, nrow = N, ncol = n)
}

optimize_beta <- function(P, L, hatS) {
  n <- ncol(hatS)
  InXL <- kronecker(diag(n), L)
  Ptemp <- pstack(P)
  term1 <- t(Ptemp) %*% InXL %*% Ptemp
  term2 <- MASS::ginv(term1) %*% t(Ptemp)
  term3 <- InXL %*% vecF(hatS)

  as.vector(term2 %*% term3)
}

GSPRegression <- function(X_noisy, P, w1, w2, max_iter = 100,
                          threshhold = 5e-4, convergence_tol = 1e-4,
                          solver = "SCS") {
  N <- nrow(X_noisy)
  n <- ncol(X_noisy)

  hatS <- X_noisy
  hatS_0 <- X_noisy
  objective <- rep(NA_real_, max_iter)
  beta <- rep(0, length(P))

  for (i in seq_len(max_iter)) {
    L <- optimize_L(hatS, P, beta, w1, w2, solver = solver)
    hatS <- optimize_hatS(L, P, beta, X_noisy, w1)
    beta <- optimize_beta(P, L, hatS)

    x <- vecF(hatS) - Pbeta(P, beta)
    D <- kronecker(diag(n), L)
    arg1 <- norm(hatS - hatS_0, type = "F")^2
    arg2 <- w1 * as.numeric(t(x) %*% D %*% x)
    arg3 <- w2 * norm(L, type = "F")^2
    objective[i] <- arg1 + arg2 + arg3

    if (i >= 3) {
      if (abs(objective[i] - objective[i - 1]) < convergence_tol) {
        objective <- objective[seq_len(i)]
        break
      }
    }
  }

  list(
    L_pruned = prune(L, threshhold),
    hatS = hatS,
    beta = beta,
    raw_L = L,
    objective = objective
  )
}

# ----------------------------------------------------------------------------
# Data prep exactly following NOAA Python logic
# ----------------------------------------------------------------------------

prepare_noaa_data <- function(filename, data_path = getwd()) {
  data <- read.csv(file.path(data_path, filename), row.names = 1, check.names = FALSE)

  signal_cols <- grep("^TMAX", names(data), value = TRUE)
  X_noisy <- as.matrix(data[, signal_cols, drop = FALSE])

  P <- list()
  regressor_prefixes <- c("ELEVATION", "TMIN")

  for (prefix in regressor_prefixes) {
    cols <- grep(paste0("^", prefix), names(data), value = TRUE)
    reg_matrix <- as.matrix(data[, cols, drop = FALSE])

    if (ncol(reg_matrix) == 1) {
      reg_matrix <- matrix(rep(reg_matrix[, 1], ncol(X_noisy)),
                           nrow = nrow(X_noisy), ncol = ncol(X_noisy))
    }

    P[[length(P) + 1]] <- reg_matrix
  }

  ones_matrix <- matrix(1, nrow = nrow(X_noisy), ncol = ncol(X_noisy))
  P[[length(P) + 1]] <- ones_matrix

  list(X_noisy = X_noisy, P = P)
}

build_partial_ground_truth <- function(dist_filename = "NOAA_dist_increased.csv",
                                       data_path = getwd()) {
  dist_matrix <- read.csv(file.path(data_path, dist_filename), row.names = 1, check.names = FALSE)
  dist_mat <- as.matrix(dist_matrix)
  N <- nrow(dist_mat)

  L_prime <- matrix(999, N, N)
  L_prime[dist_mat < 40000 & dist_mat > 0] <- -1
  L_prime[dist_mat > 193121] <- 0
  diag(L_prime) <- 0
  L_prime
}

# ----------------------------------------------------------------------------
# Run one file / all five imputations
# ----------------------------------------------------------------------------

run_one_noaa <- function(filename,
                         data_path = getwd(),
                         L_prime,
                         w1 = 0.001,
                         w2 = 0.99,
                         max_iter = 100,
                         prune_threshold = 5e-4,
                         convergence_tol = 1e-4,
                         solver = "SCS") {
  dat <- prepare_noaa_data(filename, data_path)

  fit <- GSPRegression(
    X_noisy = dat$X_noisy,
    P = dat$P,
    w1 = w1,
    w2 = w2,
    max_iter = max_iter,
    threshhold = prune_threshold,
    convergence_tol = convergence_tol,
    solver = solver
  )

  metrics <- metricsprf_partial(L_prime, fit$L_pruned)

  learned_edges <- sum(lowertrian(fit$L_pruned) != 0)
  raw_edges <- sum(lowertrian(fit$raw_L) != 0)

  list(
    file = filename,
    fit = fit,
    metrics = metrics,
    learned_edges = learned_edges,
    raw_edges = raw_edges
  )
}

run_all_noaa <- function(data_path = getwd(),
                         files = iteration_files,
                         w1 = 0.001,
                         w2 = 0.99,
                         max_iter = 100,
                         prune_threshold = 5e-4,
                         convergence_tol = 1e-4,
                         solver = "SCS") {
  L_prime <- build_partial_ground_truth("NOAA_dist_increased.csv", data_path)

  results <- lapply(files, function(f) {
    cat(sprintf("\nRunning %s\n", f))
    run_one_noaa(
      filename = f,
      data_path = data_path,
      L_prime = L_prime,
      w1 = w1,
      w2 = w2,
      max_iter = max_iter,
      prune_threshold = prune_threshold,
      convergence_tol = convergence_tol,
      solver = solver
    )
  })

  summary_df <- data.frame(
    File = sapply(results, `[[`, "file"),
    Learned_Edges = sapply(results, `[[`, "learned_edges"),
    Raw_Edges = sapply(results, `[[`, "raw_edges"),
    Precision = sapply(results, function(x) x$metrics$precision),
    Recall = sapply(results, function(x) x$metrics$recall),
    F1 = sapply(results, function(x) x$metrics$f1),
    NMI = sapply(results, function(x) x$metrics$nmi),
    Accuracy = sapply(results, function(x) x$metrics$accuracy),
    TP = sapply(results, function(x) x$metrics$tp),
    FP = sapply(results, function(x) x$metrics$fp),
    FN = sapply(results, function(x) x$metrics$fn),
    TN = sapply(results, function(x) x$metrics$tn)
  )

  pooled_df <- data.frame(
    Metric = c("Precision", "Recall", "F1", "NMI", "Accuracy"),
    Mean = c(
      mean(summary_df$Precision),
      mean(summary_df$Recall),
      mean(summary_df$F1),
      mean(summary_df$NMI),
      mean(summary_df$Accuracy)
    ),
    SD = c(
      sd(summary_df$Precision),
      sd(summary_df$Recall),
      sd(summary_df$F1),
      sd(summary_df$NMI),
      sd(summary_df$Accuracy)
    )
  )

  list(results = results, summary = summary_df, pooled = pooled_df, L_prime = L_prime)
}

# ----------------------------------------------------------------------------
# Example run
# ----------------------------------------------------------------------------

# setwd("/path/to/folder/with/NOAA/files")
# out <- run_all_noaa(data_path = getwd(), solver = "SCS")
# print(out$summary)
# print(out$pooled)
# write.csv(out$summary, "NOAA_python_matched_all_imputations.csv", row.names = FALSE)
# write.csv(out$pooled, "NOAA_python_matched_pooled.csv", row.names = FALSE)

