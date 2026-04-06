# GraphLearnR

<!-- badges: start -->
[![R-CMD-check](https://github.com/shivampawar/GraphLearnR/workflows/R-CMD-check/badge.svg)](https://github.com/shivampawar/GraphLearnR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Graph Learning with External Regressors in R**

R implementation of the GLReg (Graph Learning with External Regressors) algorithm from:

> Guo, J., Moses, S., & Wang, Z. (2023). Graph Learning From Signals With Smoothness Superimposed by Regressors. *IEEE Signal Processing Letters*, 30, 942-946.

## Overview

GraphLearnR learns graph Laplacian structure from noisy signal observations while incorporating external regressor information. This is particularly useful for:

- **Graph signal processing** with missing data
- **Network inference** from observational data
- **Spatial data analysis** with covariates
- **Time series on graphs** with external variables

### Key Features

- ✅ Learn graph structure from signals with external regressors
- ✅ Handle missing data through signal estimation
- ✅ Automatic hyperparameter optimization
- ✅ Comprehensive evaluation metrics (Precision, Recall, F-score, NMI)
- ✅ Synthetic graph generation (Barabási-Albert, Erdős-Rényi)
- ✅ Real-world applications (weather networks, sensor networks)

## Installation

### From GitHub (Recommended)

```r
# Install devtools if needed
install.packages("devtools")

# Install GraphLearnR from GitHub
devtools::install_github("shivampawar/GraphLearnR")
```

### From source

```r
# Clone the repository
git clone https://github.com/shivampawar/GraphLearnR.git

# Install
devtools::install("GraphLearnR")
```

## Quick Start

```r
library(GraphLearnR)
library(igraph)

# Generate a synthetic graph
G <- barabasi.game(20, m = 1, directed = FALSE)

# Generate signals with regressors
result <- generate_signals(G, num_signals = 100, beta = c(1, 1))

# Learn the graph structure
fit <- glreg(
  X_noisy = result$S_noisy,
  P = result$P,
  w1 = 0.001,
  w2 = 0.99
)

# Evaluate performance
L_true <- as.matrix(laplacian_matrix(G))
metrics <- evaluate_graph(L_true, fit$L)

print(metrics)
# $precision
# [1] 0.947
# 
# $recall
# [1] 0.895
# 
# $f_score
# [1] 0.920
```

## Core Functions

### `glreg()`
Main function for graph learning with external regressors.

```r
fit <- glreg(
  X_noisy,              # Observed signals (N × num_signals matrix)
  P,                    # List of regressor matrices
  w1,                   # Smoothness penalty
  w2,                   # Sparsity penalty
  threshold = 0.001,    # Edge pruning threshold
  lambda_reg = 1e-6,    # Numerical regularization
  max_iter = 100,       # Maximum iterations
  verbose = FALSE       # Print progress
)
```

**Returns:**
- `L`: Learned Laplacian matrix
- `hatS`: Estimated signals
- `beta`: Regression coefficients
- `objective`: Convergence history

### `generate_signals()`
Generate synthetic signals on a graph with external regressors.

```r
result <- generate_signals(
  G,                    # igraph object
  num_signals = 100,    # Number of signal observations
  beta = c(1, 1),       # Regression coefficients
  noise_level = 0.1,    # Noise standard deviation
  missing_rate = 0      # Proportion of missing values
)
```

### `evaluate_graph()`
Compute performance metrics between true and learned graphs.

```r
metrics <- evaluate_graph(L_true, L_learned)
```

**Returns:**
- `precision`: True positives / (True positives + False positives)
- `recall`: True positives / (True positives + False negatives)
- `f_score`: Harmonic mean of precision and recall
- `nmi`: Normalized Mutual Information

## Examples

### Example 1: Synthetic Graph (Barabási-Albert)

```r
library(GraphLearnR)
library(igraph)

# Generate BA graph
G <- barabasi.game(20, m = 1, directed = FALSE)
plot(G, main = "True Graph")

# Generate signals
result <- generate_signals(G, num_signals = 100, beta = c(1, 1, 1))

# Learn graph
fit <- glreg(result$S_noisy, result$P, w1 = 0.0003, w2 = 0.01)

# Evaluate
L_true <- as.matrix(laplacian_matrix(G))
metrics <- evaluate_graph(L_true, fit$L)

cat(sprintf("F-score: %.3f\n", metrics$f_score))
```

### Example 2: Weather Station Network

See `vignettes/weather_application.Rmd` for a complete real-world example using NOAA weather station data.

```r
# Load weather data
weather_data <- read.csv("NOAA_increased.csv", row.names = 1)

# Extract signals (TMAX) and regressors (ELEVATION, TMIN)
# ... preprocessing code ...

# Run GLReg
fit <- glreg(X_noisy, P, w1 = 0.001, w2 = 0.99)

# Evaluate against geographic groundtruth
metrics <- evaluate_graph(L_groundtruth, fit$L)
```

## Real-World Applications

This package has been successfully applied to:

1. **Weather Station Networks**: Learning temperature correlation networks from NOAA data (20 stations, 90 days)
2. **Soil Moisture Networks**: Analyzing Western Weather Group sensor data (38 stations)
3. **Synthetic Benchmarks**: BA and ER graphs with various sizes and densities

## Algorithm Details

GLReg solves the following optimization problem:

```
minimize  ||hatS - S||² + w1 · tr(hatS^T L hatS) + w2 · ||L||²_F
subject to:
  - tr(L) = N
  - L_ij ≤ 0 for i ≠ j
  - L_ii ≥ 0
  - L·1 = 0
  - hatS = P·β + smooth graph signal
```

Where:
- `S`: Observed signals
- `hatS`: Estimated signals
- `L`: Graph Laplacian
- `P`: External regressors
- `β`: Regression coefficients
- `w1, w2`: Hyperparameters

The algorithm alternates between:
1. Updating `L` (convex optimization via CVXR)
2. Updating `hatS` (closed-form solution)
3. Updating `β` (closed-form solution)

## Requirements

- R ≥ 3.5.0
- CVXR ≥ 1.8.0 (convex optimization)
- igraph ≥ 1.2.0 (graph operations)
- Matrix (sparse matrix support)

### Optional
- mice (for missing data imputation)
- ggplot2 (for visualization)

## Citation

If you use GraphLearnR in your research, please cite:

```bibtex
@mastersthesis{pawar2026graphlearnr,
  title = {GraphLearnR: An R Package for Graph Structure Learning with External Regressors},
  author = {Pawar, Shivam},
  year = {2026},
  school = {California State University, Chico},
  type = {Master's Thesis}
}

@article{guo2023graph,
  title = {Graph Learning From Signals With Smoothness Superimposed by Regressors},
  author = {Guo, Jane and Moses, Samuel and Wang, Zhihua},
  journal = {IEEE Signal Processing Letters},
  volume = {30},
  pages = {942--946},
  year = {2023}
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Dr. Jane Guo** (California State University, Chico) - Thesis advisor and GLReg algorithm developer
- **Christopher Duran** - Python implementation reference
- **Prof. Samuel Moses** - Methodology development
- **IEEE Signal Processing Society** - Original publication

## Contact

**Shivam Pawar**  
Master's Student, Data Science and Analytics  
California State University, Chico  
Email: sgpawar@csuchico.edu  
GitHub: [@shivampawar](https://github.com/shivampawar)

## Project Status

🚧 **Active Development** - Version 0.1.0 (May 2026)

This package is part of a Master's thesis project at CSU Chico.

---

**Built with ❤️ in R**
