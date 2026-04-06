#' Graph Learning with External Regressors (GLReg)
#'
#' Learns graph Laplacian from signals while accounting for external regressors.
#' Implements the alternating optimization algorithm from Guo et al. (2023).
#'
#' @param X_noisy Matrix of observed signals (N x num_signals), where N is number of nodes
#' @param P List of regressor matrices, each of dimension (N x num_signals). 
#'   The last element should be a matrix of ones for the intercept.
#' @param w1 Hyperparameter controlling smoothness penalty (higher = smoother signals)
#' @param w2 Hyperparameter controlling Laplacian sparsity (higher = sparser graph)
#' @param max_iter Maximum number of iterations (default: 100)
#' @param threshold Threshold for pruning insignificant edges (default: 0.001)
#' @param lambda_reg Regularization parameter for numerical stability (default: 1e-6).
#'   Increase if singularity errors occur with certain datasets.
#' @param tol Convergence tolerance (default: 1e-4)
#' @param verbose Logical, whether to print iteration progress (default: FALSE)
#'
#' @return A list containing:
#'   \item{L}{Learned Laplacian matrix (N x N)}
#'   \item{L_unpruned}{Unpruned Laplacian before thresholding}
#'   \item{hatS}{Estimated signals (N x num_signals)}
#'   \item{beta}{Estimated regression coefficients}
#'   \item{objective}{Vector of objective function values at each iteration}
#'   \item{iterations}{Number of iterations until convergence}
#'
#' @references
#' Guo, J., Moses, S., & Wang, Z. (2023). Graph Learning From Signals With 
#' Smoothness Superimposed by Regressors. IEEE Signal Processing Letters, 30, 942-946.
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate a random graph
#' library(igraph)
#' G <- barabasi.game(20, m = 1, directed = FALSE)
#' 
#' # Generate signals with regressors
#' result <- generate_signals(G, num_signals = 100, beta = c(1, 1))
#' S <- result$S_noisy
#' P <- result$P
#' 
#' # Learn the graph
#' fit <- glreg(S, P, w1 = 1, w2 = 0.5)
#' print(fit)
#' }
glreg <- function(X_noisy, P, w1, w2, max_iter = 100, 
                  threshold = 0.001, lambda_reg = 1e-6, tol = 1e-4, verbose = FALSE) {
  
  # Input validation
  if (!is.matrix(X_noisy)) {
    stop("X_noisy must be a matrix")
  }
  
  if (!is.list(P) || !all(sapply(P, is.matrix))) {
    stop("P must be a list of matrices")
  }
  
  if (!all(sapply(P, function(p) all(dim(p) == dim(X_noisy))))) {
    stop("All regressor matrices in P must have the same dimensions as X_noisy")
  }
  
  N <- nrow(X_noisy)
  n <- ncol(X_noisy)
  
  # Check for required packages
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required. Please install it with: install.packages('CVXR')")
  }
  
  # Initialize
  hatS <- X_noisy
  hatS_0 <- X_noisy
  beta <- matrix(0, nrow = length(P), ncol = 1)
  objective <- numeric(max_iter)
  
  if (verbose) {
    cat("Starting GLReg algorithm...\n")
    cat(sprintf("N = %d nodes, n = %d signals\n", N, n))
    cat(sprintf("w1 = %.4f, w2 = %.4f\n\n", w1, w2))
  }
  
  for (i in 1:max_iter) {
    
    # Step 1: Update Laplacian L
    L <- optimize_L(hatS, P, beta, w1, w2, N, n)
    
    # Step 2: Update signal estimate hatS
    hatS <- optimize_hatS(L, P, beta, X_noisy, w1, N, n, lambda_reg)
    
    # Step 3: Update regression coefficients beta
    beta <- optimize_beta(P, L, hatS, N, n, lambda_reg)
    
    # Calculate objective function
    arg1 <- norm(hatS - hatS_0, "F")^2
    x <- vec(hatS) - Pbeta(P, beta)
    D <- kronecker(diag(n), L)
    arg2 <- w1 * t(x) %*% D %*% x
    arg3 <- w2 * norm(L, "F")^2
    objective[i] <- arg1 + arg2[1,1] + arg3
    
    if (verbose && i %% 10 == 0) {
      cat(sprintf("Iteration %d: objective = %.6f\n", i, objective[i]))
    }
    
    # Check convergence
    if (i >= 2 && abs(objective[i] - objective[i-1]) < tol) {
      if (verbose) {
        cat(sprintf("\nConverged at iteration %d\n", i))
      }
      objective <- objective[1:i]
      break
    }
  }
  
  # Prune insignificant edges
  L_pruned <- prune_laplacian(L, threshold)
  
  # Prepare output
  result <- list(
    L = L_pruned,
    L_unpruned = L,
    hatS = hatS,
    beta = beta,
    objective = objective,
    iterations = length(objective),
    w1 = w1,
    w2 = w2,
    threshold = threshold
  )
  
  class(result) <- "glreg"
  return(result)
}


#' @keywords internal
optimize_L <- function(hatS, P, beta, w1, w2, N, n) {
  # Optimize Laplacian given hatS and beta
  # Uses CVXR for convex optimization
  # Compatible with CVXR 1.8
  
  ONE <- matrix(1, nrow = N, ncol = 1)
  ZERO <- matrix(0, nrow = N, ncol = 1)
  
  # Define optimization variable (CVXR 1.8 syntax)
  L <- CVXR::Variable(c(N, N), symmetric = TRUE)
  
  # Define constraints using CVXR 1.8 functions
  # Constraint 1: trace(L) = N
  # Constraint 2: off-diagonal elements <= 0
  # Constraint 3: diagonal elements >= 0  
  # Constraint 4: L * 1 = 0 (Laplacian property)
  
  constraints <- list(
    CVXR::matrix_trace(L) == N,           # trace(L) = N
    L %*% ONE == ZERO                     # L * 1 = 0
  )
  
  # Add elementwise constraints for diagonal and off-diagonal
  # For each element: if i == j, L[i,j] >= 0; if i != j, L[i,j] <= 0
  for (i in 1:N) {
    constraints <- c(constraints, list(L[i,i] >= 0))  # diagonal >= 0
    for (j in 1:N) {
      if (i != j) {
        constraints <- c(constraints, list(L[i,j] <= 0))  # off-diagonal <= 0
      }
    }
  }
  
  # Define objective
  x <- vec(hatS) - Pbeta(P, beta)
  
  # Reshape x back to matrix form for computation
  x_mat <- matrix(x, nrow = N, ncol = n)
  
  # Compute the smoothness term: tr(X^T L X)
  smoothness_term <- CVXR::matrix_trace(t(x_mat) %*% L %*% x_mat)
  
  objective <- CVXR::Minimize(
    w1 * smoothness_term + w2 * CVXR::sum_squares(L)
  )
  
  # Solve (CVXR 1.8 uses psolve)
  problem <- CVXR::Problem(objective, constraints)
  result <- tryCatch({
    CVXR::psolve(problem, solver = "ECOS")
  }, error = function(e) {
    # Try SCS solver as fallback
    tryCatch({
      CVXR::psolve(problem, solver = "SCS")
    }, error = function(e2) {
      stop(sprintf("Both ECOS and SCS solvers failed. Error: %s", e2$message))
    })
  })
  
  if (result$status != "optimal" && result$status != "optimal_inaccurate") {
    warning(sprintf("Laplacian optimization status: %s", result$status))
  }
  
  return(result$getValue(L))
}


#' @keywords internal
optimize_hatS <- function(L, P, beta, S, w1, N, n, lambda_reg = 1e-6) {
  # Optimize signal estimate given L and beta
  # Has closed-form solution
  # 
  # lambda_reg: Regularization parameter to ensure numerical stability
  #             Default 1e-6 works for most datasets including MICE-imputed data
  
  InXL <- kronecker(diag(n), L)
  term1 <- w1 * InXL + diag(N * n)
  
  # Add regularization to ensure invertibility (Tikhonov regularization)
  # This prevents singularity issues with different imputation methods
  term1 <- term1 + lambda_reg * diag(N * n)
  
  term2 <- solve(term1)
  term3 <- vec(S) + w1 * (InXL %*% Pbeta(P, beta))
  
  hatS_vec <- term2 %*% term3
  hatS <- matrix(hatS_vec, nrow = N, ncol = n)
  
  return(hatS)
}


#' @keywords internal
optimize_beta <- function(P, L, hatS, N, n, lambda_reg = 1e-6) {
  # Optimize regression coefficients given L and hatS
  # Has closed-form solution (weighted least squares)
  
  InXL <- kronecker(diag(n), L)
  Pstack <- pstack(P)
  
  term1 <- t(Pstack) %*% InXL %*% Pstack
  
  # Add regularization to ensure invertibility
  term1 <- term1 + lambda_reg * diag(ncol(Pstack))
  
  term2 <- solve(term1) %*% t(Pstack)
  term3 <- InXL %*% vec(hatS)
  
  beta <- term2 %*% term3
  
  return(beta)
}


#' @keywords internal
vec <- function(S) {
  # Vectorize matrix column-wise (Fortran order)
  return(as.vector(S))
}


#' @keywords internal
Pbeta <- function(P, beta) {
  # Compute P*beta where P is list of matrices and beta is coefficient vector
  R <- Reduce(`+`, mapply(function(p, b) p * b, P, as.vector(beta), SIMPLIFY = FALSE))
  return(vec(R))
}


#' @keywords internal
pstack <- function(P) {
  # Stack list of regressor matrices into a single matrix
  # Each matrix is vectorized column-wise and stacked as columns
  Pvec <- lapply(P, vec)
  Pstack <- do.call(cbind, Pvec)
  return(Pstack)
}


#' @keywords internal
prune_laplacian <- function(L, threshold) {
  # Set edge weights below threshold to zero
  L_pruned <- L
  L_pruned[abs(L_pruned) < threshold] <- 0
  return(L_pruned)
}


#' Print method for glreg objects
#' @param x A glreg object
#' @param ... Additional arguments (not used)
#' @export
print.glreg <- function(x, ...) {
  cat("GLReg: Graph Learning with External Regressors\n")
  cat("==============================================\n\n")
  cat(sprintf("Number of nodes: %d\n", nrow(x$L)))
  cat(sprintf("Number of signals: %d\n", ncol(x$hatS)))
  cat(sprintf("Number of regressors: %d\n", length(x$beta)))
  cat(sprintf("Iterations: %d\n", x$iterations))
  cat(sprintf("Final objective: %.6f\n\n", tail(x$objective, 1)))
  cat(sprintf("Hyperparameters: w1 = %.4f, w2 = %.4f\n", x$w1, x$w2))
  cat(sprintf("Edge threshold: %.4f\n\n", x$threshold))
  
  # Count edges
  n_edges <- sum(x$L[lower.tri(x$L)] != 0)
  n_possible <- choose(nrow(x$L), 2)
  cat(sprintf("Learned edges: %d / %d (%.1f%% density)\n", 
              n_edges, n_possible, 100 * n_edges / n_possible))
  
  cat("\nRegression coefficients:\n")
  print(x$beta)
}


#' Summary method for glreg objects
#' @param object A glreg object
#' @param ... Additional arguments (not used)
#' @export
summary.glreg <- function(object, ...) {
  cat("GLReg Summary\n")
  cat("=============\n\n")
  print(object)
  
  cat("\n\nObjective function convergence:\n")
  cat(sprintf("Initial: %.6f\n", object$objective[1]))
  cat(sprintf("Final:   %.6f\n", tail(object$objective, 1)))
  cat(sprintf("Change:  %.6f\n", object$objective[1] - tail(object$objective, 1)))
  
  invisible(object)
}
