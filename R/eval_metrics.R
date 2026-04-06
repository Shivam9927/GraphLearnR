#' Evaluate Learned Graph Against Ground Truth
#'
#' Computes precision, recall, F-score, and normalized mutual information (NMI)
#' for comparing a learned graph Laplacian against ground truth.
#'
#' @param L_true True Laplacian matrix (N x N)
#' @param L_learned Learned Laplacian matrix (N x N)
#' @param use_lower_tri Logical, whether to use only lower triangular part (default: TRUE)
#'
#' @return A list containing:
#'   \item{precision}{Precision = TP / (TP + FP)}
#'   \item{recall}{Recall = TP / (TP + FN)}
#'   \item{f_score}{F-score = 2 * precision * recall / (precision + recall)}
#'   \item{nmi}{Normalized Mutual Information}
#'   \item{TP}{True positives (correct edges)}
#'   \item{FP}{False positives (incorrect edges)}
#'   \item{FN}{False negatives (missed edges)}
#'   \item{TN}{True negatives (correct non-edges)}
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate a graph and signals
#' G <- sample_pa(20, m = 1, directed = FALSE)
#' result <- generate_signals(G, 100, beta = c(1, 0.5))
#' 
#' # Learn the graph
#' fit <- glreg(result$S_noisy, result$P, w1 = 1, w2 = 0.5)
#' 
#' # Evaluate
#' metrics <- evaluate_graph(result$L, fit$L)
#' print(metrics)
#' }
evaluate_graph <- function(L_true, L_learned, use_lower_tri = TRUE) {
  
  # Input validation
  if (!is.matrix(L_true) || !is.matrix(L_learned)) {
    stop("Both L_true and L_learned must be matrices")
  }
  
  if (!all(dim(L_true) == dim(L_learned))) {
    stop("L_true and L_learned must have the same dimensions")
  }
  
  N <- nrow(L_true)
  
  if (use_lower_tri) {
    # Extract lower triangular part (excluding diagonal)
    lower_tri_idx <- lower.tri(L_true, diag = FALSE)
    edges_true <- L_true[lower_tri_idx] != 0
    edges_learned <- L_learned[lower_tri_idx] != 0
  } else {
    # Use all off-diagonal elements
    diag_mask <- diag(N) == 1
    edges_true <- (L_true != 0) & !diag_mask
    edges_learned <- (L_learned != 0) & !diag_mask
    edges_true <- as.vector(edges_true)
    edges_learned <- as.vector(edges_learned)
  }
  
  # Calculate confusion matrix
  TP <- sum(edges_true & edges_learned)    # True positives
  FP <- sum(!edges_true & edges_learned)   # False positives
  FN <- sum(edges_true & !edges_learned)   # False negatives
  TN <- sum(!edges_true & !edges_learned)  # True negatives
  
  # Calculate metrics
  if ((TP + FP) > 0) {
    precision <- TP / (TP + FP)
  } else {
    precision <- 0
    warning("No edges learned (TP + FP = 0), precision set to 0")
  }
  
  if ((TP + FN) > 0) {
    recall <- TP / (TP + FN)
  } else {
    recall <- 0
    warning("No true edges (TP + FN = 0), recall set to 0")
  }
  
  if ((precision + recall) > 0) {
    f_score <- 2 * precision * recall / (precision + recall)
  } else {
    f_score <- 0
  }
  
  # Calculate NMI (Normalized Mutual Information)
  nmi <- calculate_nmi(edges_true, edges_learned)
  
  # Create result object
  result <- list(
    precision = precision,
    recall = recall,
    f_score = f_score,
    nmi = nmi,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    n_true_edges = sum(edges_true),
    n_learned_edges = sum(edges_learned),
    n_possible_edges = length(edges_true)
  )
  
  class(result) <- "graph_evaluation"
  return(result)
}


#' Calculate Normalized Mutual Information
#'
#' Computes NMI between true and predicted edge indicators
#'
#' @param y_true Logical vector of true edge indicators
#' @param y_pred Logical vector of predicted edge indicators
#'
#' @return Numeric NMI value between 0 and 1
#' @keywords internal
calculate_nmi <- function(y_true, y_pred) {
  
  # Convert to numeric
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  
  # Calculate contingency table
  n00 <- sum(y_true == 0 & y_pred == 0)
  n01 <- sum(y_true == 0 & y_pred == 1)
  n10 <- sum(y_true == 1 & y_pred == 0)
  n11 <- sum(y_true == 1 & y_pred == 1)
  
  n <- length(y_true)
  
  # Calculate marginals
  n0_true <- n00 + n01
  n1_true <- n10 + n11
  n0_pred <- n00 + n10
  n1_pred <- n01 + n11
  
  # Calculate entropies
  # H(Y_true)
  p0_true <- n0_true / n
  p1_true <- n1_true / n
  h_true <- 0
  if (p0_true > 0) h_true <- h_true - p0_true * log2(p0_true)
  if (p1_true > 0) h_true <- h_true - p1_true * log2(p1_true)
  
  # H(Y_pred)
  p0_pred <- n0_pred / n
  p1_pred <- n1_pred / n
  h_pred <- 0
  if (p0_pred > 0) h_pred <- h_pred - p0_pred * log2(p0_pred)
  if (p1_pred > 0) h_pred <- h_pred - p1_pred * log2(p1_pred)
  
  # H(Y_true, Y_pred) - joint entropy
  h_joint <- 0
  if (n00 > 0) h_joint <- h_joint - (n00/n) * log2(n00/n)
  if (n01 > 0) h_joint <- h_joint - (n01/n) * log2(n01/n)
  if (n10 > 0) h_joint <- h_joint - (n10/n) * log2(n10/n)
  if (n11 > 0) h_joint <- h_joint - (n11/n) * log2(n11/n)
  
  # Mutual information: I(Y_true; Y_pred) = H(Y_true) + H(Y_pred) - H(Y_true, Y_pred)
  mi <- h_true + h_pred - h_joint
  
  # Normalized MI: NMI = 2*I / (H(Y_true) + H(Y_pred))
  if ((h_true + h_pred) > 0) {
    nmi <- 2 * mi / (h_true + h_pred)
  } else {
    nmi <- 0
  }
  
  # Ensure NMI is in [0, 1]
  nmi <- max(0, min(1, nmi))
  
  return(nmi)
}


#' Print method for graph_evaluation objects
#' @param x A graph_evaluation object
#' @param ... Additional arguments (not used)
#' @export
print.graph_evaluation <- function(x, ...) {
  cat("Graph Learning Evaluation Metrics\n")
  cat("==================================\n\n")
  
  cat("Edge Counts:\n")
  cat(sprintf("  True edges:    %d\n", x$n_true_edges))
  cat(sprintf("  Learned edges: %d\n", x$n_learned_edges))
  cat(sprintf("  Possible edges: %d\n\n", x$n_possible_edges))
  
  cat("Confusion Matrix:\n")
  cat(sprintf("  True Positives (TP):  %d\n", x$TP))
  cat(sprintf("  False Positives (FP): %d\n", x$FP))
  cat(sprintf("  False Negatives (FN): %d\n", x$FN))
  cat(sprintf("  True Negatives (TN):  %d\n\n", x$TN))
  
  cat("Performance Metrics:\n")
  cat(sprintf("  Precision: %.4f\n", x$precision))
  cat(sprintf("  Recall:    %.4f\n", x$recall))
  cat(sprintf("  F-score:   %.4f\n", x$f_score))
  cat(sprintf("  NMI:       %.4f\n", x$nmi))
  
  if (x$precision == 1.0 && x$recall < 0.5) {
    cat("\nNote: High precision but low recall suggests conservative learning.\n")
    cat("Consider reducing sparsity penalty (w2) or threshold.\n")
  }
  
  if (x$recall == 1.0 && x$precision < 0.5) {
    cat("\nNote: High recall but low precision suggests over-fitting.\n")
    cat("Consider increasing sparsity penalty (w2) or threshold.\n")
  }
}


#' Compare multiple graph learning results
#'
#' Evaluates multiple learned graphs against ground truth and compares them
#'
#' @param L_true True Laplacian matrix
#' @param L_list List of learned Laplacian matrices to compare
#' @param method_names Character vector of method names (optional)
#'
#' @return A data frame with comparison metrics
#' @export
compare_methods <- function(L_true, L_list, method_names = NULL) {
  
  if (!is.list(L_list)) {
    stop("L_list must be a list of matrices")
  }
  
  if (is.null(method_names)) {
    method_names <- paste0("Method_", 1:length(L_list))
  }
  
  if (length(method_names) != length(L_list)) {
    stop("Length of method_names must match length of L_list")
  }
  
  # Evaluate each method
  results <- lapply(L_list, function(L) evaluate_graph(L_true, L))
  
  # Create comparison data frame
  comparison <- data.frame(
    Method = method_names,
    Precision = sapply(results, function(r) r$precision),
    Recall = sapply(results, function(r) r$recall),
    F_score = sapply(results, function(r) r$f_score),
    NMI = sapply(results, function(r) r$nmi),
    Learned_Edges = sapply(results, function(r) r$n_learned_edges),
    True_Edges = sapply(results, function(r) r$n_true_edges)
  )
  
  # Sort by F-score (descending)
  comparison <- comparison[order(-comparison$F_score), ]
  
  class(comparison) <- c("method_comparison", "data.frame")
  return(comparison)
}


#' Print method for method_comparison objects
#' @param x A method_comparison object
#' @param ... Additional arguments (not used)
#' @export
print.method_comparison <- function(x, ...) {
  cat("Method Comparison Results\n")
  cat("=========================\n\n")
  
  # Print as formatted table
  print(x, row.names = FALSE, digits = 4)
  
  cat("\n")
  best_method <- x$Method[1]
  best_f_score <- x$F_score[1]
  cat(sprintf("Best method: %s (F-score = %.4f)\n", best_method, best_f_score))
}


#' Plot comparison of methods
#' @param x A method_comparison object
#' @param ... Additional arguments (not used)
#' @export
plot.method_comparison <- function(x, ...) {
  par(mfrow = c(1, 3))
  
  # Plot 1: Precision
  barplot(x$Precision, names.arg = x$Method,
          main = "Precision", ylab = "Precision",
          col = "steelblue", las = 2, ylim = c(0, 1))
  abline(h = 0.5, lty = 2, col = "red")
  
  # Plot 2: Recall
  barplot(x$Recall, names.arg = x$Method,
          main = "Recall", ylab = "Recall",
          col = "forestgreen", las = 2, ylim = c(0, 1))
  abline(h = 0.5, lty = 2, col = "red")
  
  # Plot 3: F-score
  barplot(x$F_score, names.arg = x$Method,
          main = "F-score", ylab = "F-score",
          col = "darkorange", las = 2, ylim = c(0, 1))
  abline(h = 0.5, lty = 2, col = "red")
  
  par(mfrow = c(1, 1))
}


#' Calculate metrics for a single glreg result object
#'
#' Convenience function to evaluate a glreg object directly
#'
#' @param glreg_result A glreg object (output from glreg())
#' @param L_true True Laplacian matrix
#'
#' @return A graph_evaluation object
#' @export
evaluate_glreg <- function(glreg_result, L_true) {
  if (!inherits(glreg_result, "glreg")) {
    stop("glreg_result must be a glreg object")
  }
  
  evaluate_graph(L_true, glreg_result$L)
}
