suppressPackageStartupMessages({
  library(CVXR)
  library(MASS)
})

# =============================================================================
# NOAA Graph Learning in R
# Rewritten to follow the methodology described in:
# "Graph Signal Prediction using Learned Laplacian with Application to Weather Data"
#
# Paper-aligned choices for NOAA:
# - Signal: TMAX
# - Regressors: ELEVATION, TMIN
# - 5 imputed datasets
# - Partial ground-truth graph built from station distances
# - Hyperparameter selection by grid search on F1
# - Final pruning threshold for NOAA: 0.005
#
# Note:
# This script follows the paper methodology, but exact numerical agreement with
# the Python implementation may still vary because CVXR/CVXPY solvers can differ.
# =============================================================================

iteration_files <- c(
  "NOAA_increased.csv",
  "NOAA_iteration2.csv",
  "NOAA_iteration3.csv",
  "NOAA_iteration4.csv",
  "NOAA_iteration5.csv"
)

default_prune_threshold <- 0.005
default_gamma1_grid <- c(1e-4, 5e-4, 9e-4, 1e-3, 2e-3, 5e-3)
default_gamma2_grid <- c(0.50, 0.80, 0.90, 0.95, 0.99, 1.05, 1.15)

lowertrian <- function(Lap) {
  Lap[lower.tri(Lap, diag = FALSE)]
}

vecF <- function(S) {
  as.vector(S)
}

pstack <- function(P) {
  do.call(cbind, lapply(P, vecF))
}

Pbeta <- function(P, beta) {
  as.vector(pstack(P) %*% as.numeric(beta))
}

prune <- function(L, threshold) {
  out <- L
  out[abs(out) < threshold] <- 0
  out
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
  Lprime_tmp <- lowertrian(Lprime)
  L_tmp <- lowertrian(L)

  known_idx <- which(Lprime_tmp != 999)
  if (length(known_idx) == 0) {
    return(list(
      precision = NA_real_, recall = NA_real_, f1 = NA_real_,
      nmi = NA_real_, accuracy = NA_real_,
      tp = NA_integer_, fp = NA_integer_, fn = NA_integer_, tn = NA_integer_
    ))
  }

  true_edges <- (Lprime_tmp[known_idx] != 0)
  learned_edges <- (L_tmp[known_idx] != 0)

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

build_partial_ground_truth <- function(dist_filename = "NOAA_dist_increased.csv",
                                       data_path = getwd()) {
  dist_matrix <- read.csv(file.path(data_path, dist_filename),
                          row.names = 1, check.names = FALSE)
  dist_mat <- as.matrix(dist_matrix)

  L_prime <- matrix(999, nrow(dist_mat), ncol(dist_mat))
  L_prime[dist_mat < 40000 & dist_mat > 0] <- -1
  L_prime[dist_mat > 193121] <- 0
  diag(L_prime) <- 0
  L_prime
}

prepare_noaa_data <- function(filename, data_path = getwd()) {
  data <- read.csv(file.path(data_path, filename),
                   row.names = 1, check.names = FALSE)

  tmax_cols <- grep("^TMAX", names(data), value = TRUE)
  tmin_cols <- grep("^TMIN", names(data), value = TRUE)

  X_noisy <- as.matrix(data[, tmax_cols, drop = FALSE])
  TMIN <- as.matrix(data[, tmin_cols, drop = FALSE])
  elevation <- data$ELEVATION

  N <- nrow(X_noisy)
  n <- ncol(X_noisy)

  P1 <- matrix(rep(elevation, n), nrow = N, ncol = n)
  P2 <- TMIN
  P3 <- matrix(1, nrow = N, ncol = n)

  list(
    X_noisy = X_noisy,
    P = list(P1, P2, P3)
  )
}

optimize_L <- function(hatS, P, beta, gamma1, gamma2, solver = "SCS") {
  N <- nrow(hatS)
  L <- Variable(c(N, N), symmetric = TRUE)

  R <- hatS - P[[1]] * beta[1] - P[[2]] * beta[2] - P[[3]] * beta[3]
  S_mat <- R %*% t(R)
  one_vec <- matrix(1, nrow = N, ncol = 1)

  constraints <- list(
    sum(diag(L)) == N,
    L %*% one_vec == 0,
    diag(L) >= 0
  )

  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (i != j) {
        constraints <- c(constraints, L[i, j] <= 0)
      }
    }
  }

  L_vec <- vec(L)
  linear_term <- sum_entries(multiply(L, S_mat))
  quad_term <- sum_squares(L_vec)

  objective <- Minimize(gamma1 * linear_term + gamma2 * quad_term)

  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = solver, verbose = FALSE)

  if (result$status %in% c("infeasible", "unbounded")) {
    stop(sprintf("optimize_L failed with solver status: %s", result$status))
  }

  out <- result$getValue(L)
  out <- (out + t(out)) / 2
  out
}

optimize_hats <- function(L, P, beta, S, gamma1) {
  N <- nrow(S)
  n <- ncol(S)

  A <- diag(N * n) + gamma1 * kronecker(diag(n), L)
  b <- vecF(S) + gamma1 * (kronecker(diag(n), L) %*% Pbeta(P, beta))

  hats_vec <- MASS::ginv(A) %*% b
  matrix(hats_vec, nrow = N, ncol = n)
}

optimize_beta <- function(P, L, hats) {
  n <- ncol(hats)
  InXL <- kronecker(diag(n), L)
  Ptemp <- pstack(P)

  A <- t(Ptemp) %*% InXL %*% Ptemp
  b <- t(Ptemp) %*% InXL %*% vecF(hats)

  as.vector(MASS::ginv(A) %*% b)
}

GSPRegression <- function(X_noisy, P, gamma1, gamma2,
                          max_iter = 100,
                          prune_threshold = default_prune_threshold,
                          convergence_tol = 1e-4,
                          solver = "SCS") {
  N <- nrow(X_noisy)
  n <- ncol(X_noisy)

  hats <- X_noisy
  hats0 <- X_noisy
  beta <- rep(0, length(P))
  objective <- rep(NA_real_, max_iter)
  raw_L <- matrix(0, N, N)

  for (i in seq_len(max_iter)) {
    raw_L <- optimize_L(hats, P, beta, gamma1, gamma2, solver = solver)
    hats <- optimize_hats(raw_L, P, beta, X_noisy, gamma1)
    beta <- optimize_beta(P, raw_L, hats)

    x <- vecF(hats) - Pbeta(P, beta)
    D <- kronecker(diag(n), raw_L)
    objective[i] <- norm(hats - hats0, type = "F")^2 +
      gamma1 * as.numeric(t(x) %*% D %*% x) +
      gamma2 * norm(raw_L, type = "F")^2

    if (i >= 3 && abs(objective[i] - objective[i - 1]) < convergence_tol) {
      objective <- objective[seq_len(i)]
      break
    }
  }

  list(
    L_pruned = prune(raw_L, prune_threshold),
    hatS = hats,
    beta = beta,
    raw_L = raw_L,
    objective = objective
  )
}

score_one_setting <- function(filename, data_path, L_prime,
                              gamma1, gamma2,
                              max_iter, prune_threshold,
                              convergence_tol, solver) {
  dat <- prepare_noaa_data(filename, data_path)

  fit <- GSPRegression(
    X_noisy = dat$X_noisy,
    P = dat$P,
    gamma1 = gamma1,
    gamma2 = gamma2,
    max_iter = max_iter,
    prune_threshold = prune_threshold,
    convergence_tol = convergence_tol,
    solver = solver
  )

  metrics <- metricsprf_partial(L_prime, fit$L_pruned)

  data.frame(
    File = filename,
    gamma1 = gamma1,
    gamma2 = gamma2,
    Precision = metrics$precision,
    Recall = metrics$recall,
    F1 = metrics$f1,
    NMI = metrics$nmi,
    Accuracy = metrics$accuracy,
    TP = metrics$tp,
    FP = metrics$fp,
    FN = metrics$fn,
    TN = metrics$tn,
    Learned_Edges = sum(lowertrian(fit$L_pruned) != 0),
    Raw_Edges = sum(lowertrian(fit$raw_L) != 0),
    Iterations = length(fit$objective)
  )
}

grid_search_one_file <- function(filename,
                                 data_path = getwd(),
                                 L_prime,
                                 gamma1_grid = default_gamma1_grid,
                                 gamma2_grid = default_gamma2_grid,
                                 max_iter = 100,
                                 prune_threshold = default_prune_threshold,
                                 convergence_tol = 1e-4,
                                 solver = "SCS") {
  results <- list()
  k <- 1L

  for (g1 in gamma1_grid) {
    for (g2 in gamma2_grid) {
      cat(sprintf("  trying gamma1=%.6f gamma2=%.6f\n", g1, g2))
      results[[k]] <- score_one_setting(
        filename = filename,
        data_path = data_path,
        L_prime = L_prime,
        gamma1 = g1,
        gamma2 = g2,
        max_iter = max_iter,
        prune_threshold = prune_threshold,
        convergence_tol = convergence_tol,
        solver = solver
      )
      k <- k + 1L
    }
  }

  grid_df <- do.call(rbind, results)
  grid_df <- grid_df[order(-grid_df$F1, -grid_df$Recall, -grid_df$Precision), ]
  rownames(grid_df) <- NULL

  list(best = grid_df[1, , drop = FALSE], full_grid = grid_df)
}

run_noaa_paper_method <- function(data_path = getwd(),
                                  files = iteration_files,
                                  gamma1_grid = default_gamma1_grid,
                                  gamma2_grid = default_gamma2_grid,
                                  max_iter = 100,
                                  prune_threshold = default_prune_threshold,
                                  convergence_tol = 1e-4,
                                  solver = "SCS") {
  L_prime <- build_partial_ground_truth(data_path = data_path)

  best_rows <- list()
  full_grids <- list()

  for (i in seq_along(files)) {
    cat(sprintf("\nGrid search for %s\n", files[i]))
    out_i <- grid_search_one_file(
      filename = files[i],
      data_path = data_path,
      L_prime = L_prime,
      gamma1_grid = gamma1_grid,
      gamma2_grid = gamma2_grid,
      max_iter = max_iter,
      prune_threshold = prune_threshold,
      convergence_tol = convergence_tol,
      solver = solver
    )
    best_rows[[i]] <- out_i$best
    full_grids[[i]] <- out_i$full_grid
  }

  summary_df <- do.call(rbind, best_rows)
  rownames(summary_df) <- NULL

  pooled <- data.frame(
    Metric = c("Precision", "Recall", "F1", "NMI", "Accuracy"),
    Mean = c(
      mean(summary_df$Precision, na.rm = TRUE),
      mean(summary_df$Recall, na.rm = TRUE),
      mean(summary_df$F1, na.rm = TRUE),
      mean(summary_df$NMI, na.rm = TRUE),
      mean(summary_df$Accuracy, na.rm = TRUE)
    ),
    SD = c(
      sd(summary_df$Precision, na.rm = TRUE),
      sd(summary_df$Recall, na.rm = TRUE),
      sd(summary_df$F1, na.rm = TRUE),
      sd(summary_df$NMI, na.rm = TRUE),
      sd(summary_df$Accuracy, na.rm = TRUE)
    )
  )

  selected_hyperparameters <- data.frame(
    File = summary_df$File,
    gamma1 = summary_df$gamma1,
    gamma2 = summary_df$gamma2
  )

  list(
    summary = summary_df,
    pooled = pooled,
    selected_hyperparameters = selected_hyperparameters,
    full_grids = full_grids,
    L_prime = L_prime
  )
}
