#' Generate Signals from Graph with External Regressors
#'
#' Generates graph signals according to the model: S = Uf + P*beta + noise
#' where U are eigenvectors of the graph Laplacian, f are latent factors,
#' P are external regressors, and beta are regression coefficients.
#'
#' @param graph An igraph object representing the graph structure
#' @param num_signals Integer, number of signal observations to generate
#' @param beta Numeric vector of regression coefficients (including intercept)
#' @param P Optional list of regressor matrices. If NULL, will generate random normal regressors.
#'   Each matrix should be (num_nodes x num_signals). The first should be a matrix of 1's for intercept.
#' @param mu Mean for generating latent factors (default: 0)
#' @param sigma Standard deviation for observation noise (default: 0.5)
#' @param regressor_mean Mean for generating random regressors (default: 0)
#' @param regressor_sd Standard deviation for generating random regressors (default: 10)
#' @param seed Random seed for reproducibility
#'
#' @return A list containing:
#'   \item{S_clean}{Clean signals without noise (N x num_signals)}
#'   \item{S_noisy}{Noisy signals (N x num_signals)}
#'   \item{P}{List of regressor matrices used}
#'   \item{L}{Normalized Laplacian matrix used}
#'   \item{beta}{Regression coefficients used}
#'   \item{U}{Eigenvectors of Laplacian}
#'   \item{Lambda}{Eigenvalues of Laplacian}
#'
#' @export
#' @examples
#' \dontrun{
#' library(igraph)
#' G <- sample_pa(20, m = 1, directed = FALSE)
#' result <- generate_signals(G, num_signals = 100, beta = c(1, 0.5))
#' plot(result$S_noisy[,1], main = "First signal observation")
#' }
generate_signals <- function(graph, 
                             num_signals, 
                             beta, 
                             P = NULL,
                             mu = 0, 
                             sigma = 0.5,
                             regressor_mean = 0,
                             regressor_sd = 10,
                             seed = NULL) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it.")
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Get the graph Laplacian
  L <- as.matrix(igraph::laplacian_matrix(graph))
  
  # Number of vertices
  N <- nrow(L)
  
  # Normalize Laplacian: trace(L) = N
  NormL <- (N / sum(diag(L))) * L
  
  # Get eigendecomposition
  eigen_decomp <- eigen(NormL)
  U <- eigen_decomp$vectors  # Eigenvectors
  Lambda <- eigen_decomp$values  # Eigenvalues
  
  # Compute covariance matrix for latent factors (pseudoinverse of eigenvalues)
  # Use pinv to handle zero/small eigenvalues
  Lambda_inv <- ifelse(abs(Lambda) > 1e-10, 1/Lambda, 0)
  cov_h <- diag(Lambda_inv)
  
  # Generate latent factors from multivariate normal
  my_mean <- rep(mu, N)
  
  # Generate smooth component using eigenvector decomposition
  # h ~ N(my_mean, cov_h)
  h_samples <- matrix(0, nrow = N, ncol = num_signals)
  
  # Generate latent factors with correct covariance structure
  # We need h such that Cov(h) = U * diag(Lambda_inv) * U^T
  for (i in 1:num_signals) {
    # Generate standard normal
    z <- rnorm(N)
    # Transform: h = U * sqrt(Lambda_inv) * z + my_mean
    sqrt_cov <- sqrt(pmax(Lambda_inv, 0))  # Ensure non-negative
    h_samples[, i] <- my_mean + U %*% (sqrt_cov * z)
  }
  
  # Graph smooth component
  X_smooth <- h_samples
  
  # Generate or use provided regressors
  if (is.null(P)) {
    # Generate default regressors
    # First regressor: intercept (all ones)
    X1 <- matrix(1, nrow = N, ncol = num_signals)
    
    # Additional regressors: random normal
    P <- list(X1)
    
    for (j in 2:length(beta)) {
      Xj <- matrix(rnorm(N * num_signals, 
                         mean = regressor_mean, 
                         sd = regressor_sd), 
                   nrow = N, ncol = num_signals)
      P[[j]] <- Xj
    }
  }
  
  # Validate P and beta dimensions
  if (length(P) != length(beta)) {
    stop(sprintf("Number of regressors (%d) must match length of beta (%d)", 
                 length(P), length(beta)))
  }
  
  # Compute regressor component: R = sum(beta[i] * P[[i]])
  R <- matrix(0, nrow = N, ncol = num_signals)
  for (i in 1:length(beta)) {
    R <- R + beta[i] * P[[i]]
  }
  
  # Generate clean signals: S = Uf + R
  S_clean <- X_smooth + R
  
  # Add observation noise
  noise <- matrix(rnorm(N * num_signals, mean = 0, sd = sigma), 
                  nrow = N, ncol = num_signals)
  S_noisy <- S_clean + noise
  
  # Return results
  result <- list(
    S_clean = S_clean,
    S_noisy = S_noisy,
    P = P,
    L = NormL,
    beta = beta,
    U = U,
    Lambda = Lambda,
    graph = graph,
    parameters = list(
      num_signals = num_signals,
      mu = mu,
      sigma = sigma,
      seed = seed
    )
  )
  
  class(result) <- "glreg_signals"
  return(result)
}


#' Print method for glreg_signals objects
#' @param x A glreg_signals object
#' @param ... Additional arguments (not used)
#' @export
print.glreg_signals <- function(x, ...) {
  cat("Graph Signal Generation Results\n")
  cat("================================\n\n")
  cat(sprintf("Graph: %d nodes, %d edges\n", 
              nrow(x$S_clean), 
              igraph::ecount(x$graph)))
  cat(sprintf("Signals: %d observations\n", ncol(x$S_clean)))
  cat(sprintf("Regressors: %d (including intercept)\n", length(x$P)))
  cat(sprintf("Noise SD: %.3f\n\n", x$parameters$sigma))
  
  cat("Regression coefficients (beta):\n")
  for (i in 1:length(x$beta)) {
    cat(sprintf("  beta[%d] = %.4f\n", i, x$beta[i]))
  }
  
  cat("\nSignal statistics:\n")
  cat(sprintf("  Clean signal range: [%.3f, %.3f]\n", 
              min(x$S_clean), max(x$S_clean)))
  cat(sprintf("  Noisy signal range: [%.3f, %.3f]\n", 
              min(x$S_noisy), max(x$S_noisy)))
  cat(sprintf("  SNR (approx): %.2f dB\n", 
              10 * log10(var(as.vector(x$S_clean)) / 
                        var(as.vector(x$S_noisy - x$S_clean)))))
}


#' Plot method for glreg_signals objects
#' @param x A glreg_signals object
#' @param ... Additional arguments (not used)
#' @export
plot.glreg_signals <- function(x, ...) {
  par(mfrow = c(2, 2))
  
  # Plot 1: First clean signal
  plot(x$S_clean[, 1], 
       type = "b", pch = 19, col = "blue",
       main = "First Clean Signal",
       xlab = "Node", ylab = "Signal Value",
       las = 1)
  grid()
  
  # Plot 2: First noisy signal
  plot(x$S_noisy[, 1], 
       type = "b", pch = 19, col = "red",
       main = "First Noisy Signal",
       xlab = "Node", ylab = "Signal Value",
       las = 1)
  grid()
  
  # Plot 3: Signal heatmap (clean)
  N <- nrow(x$S_clean)
  n <- min(50, ncol(x$S_clean))  # Show first 50 signals max
  image(1:n, 1:N, t(x$S_clean[, 1:n]),
        main = "Clean Signals Heatmap",
        xlab = "Signal Index", ylab = "Node",
        col = heat.colors(20),
        las = 1)
  
  # Plot 4: Graph structure
  if (requireNamespace("igraph", quietly = TRUE)) {
    igraph::plot.igraph(x$graph,
                       vertex.size = 15,
                       vertex.color = "lightblue",
                       vertex.label.cex = 0.8,
                       edge.width = 2,
                       main = "Graph Structure",
                       layout = igraph::layout_with_fr)
  }
  
  par(mfrow = c(1, 1))
}


#' Generate signals with specific regressor distributions
#'
#' Helper function to generate signals with different regressor types
#' (normal, binomial, exponential, etc.)
#'
#' @param graph An igraph object
#' @param num_signals Number of signals to generate
#' @param beta Regression coefficients
#' @param regressor_type Character, one of: "normal", "binomial", "exponential", 
#'   "poisson", "gamma", "uniform"
#' @param mu Mean for latent factors (default: 0)
#' @param sigma Noise standard deviation (default: 0.5)
#' @param seed Random seed
#' @param ... Additional parameters for the specific distribution
#'
#' @return A glreg_signals object
#' @export
generate_signals_with_distribution <- function(graph, 
                                              num_signals,
                                              beta,
                                              regressor_type = "normal",
                                              mu = 0,
                                              sigma = 0.5,
                                              seed = NULL,
                                              ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  N <- igraph::vcount(graph)
  
  # Intercept (all ones)
  X1 <- matrix(1, nrow = N, ncol = num_signals)
  
  # Generate regressor based on type
  X2 <- switch(regressor_type,
    "normal" = {
      args <- list(...)
      mean <- ifelse(is.null(args$mean), 0, args$mean)
      sd <- ifelse(is.null(args$sd), 10, args$sd)
      matrix(rnorm(N * num_signals, mean = mean, sd = sd), 
             nrow = N, ncol = num_signals)
    },
    "binomial" = {
      args <- list(...)
      prob <- ifelse(is.null(args$prob), 0.5, args$prob)
      matrix(rbinom(N * num_signals, size = 1, prob = prob),
             nrow = N, ncol = num_signals)
    },
    "exponential" = {
      args <- list(...)
      rate <- ifelse(is.null(args$rate), 1, args$rate)
      matrix(rexp(N * num_signals, rate = rate),
             nrow = N, ncol = num_signals)
    },
    "poisson" = {
      args <- list(...)
      lambda <- ifelse(is.null(args$lambda), 1, args$lambda)
      matrix(rpois(N * num_signals, lambda = lambda),
             nrow = N, ncol = num_signals)
    },
    "gamma" = {
      args <- list(...)
      shape <- ifelse(is.null(args$shape), 1, args$shape)
      scale <- ifelse(is.null(args$scale), 1, args$scale)
      matrix(rgamma(N * num_signals, shape = shape, scale = scale),
             nrow = N, ncol = num_signals)
    },
    "uniform" = {
      args <- list(...)
      min_val <- ifelse(is.null(args$min), 0, args$min)
      max_val <- ifelse(is.null(args$max), 1, args$max)
      matrix(runif(N * num_signals, min = min_val, max = max_val),
             nrow = N, ncol = num_signals)
    },
    stop(sprintf("Unknown regressor_type: %s", regressor_type))
  )
  
  P <- list(X1, X2)
  
  # Generate signals using main function
  generate_signals(graph, num_signals, beta, P = P, 
                  mu = mu, sigma = sigma, seed = seed)
}
