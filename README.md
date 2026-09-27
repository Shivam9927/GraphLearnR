# GraphLearnR

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Graph Learning with External Regressors in R**

GraphLearnR is an R implementation of the GLReg (Graph Learning with External Regressors) algorithm from:

> Guo, J., Moses, S., & Wang, Z. (2023). Graph Learning From Signals With Smoothness Superimposed by Regressors. *IEEE Signal Processing Letters*, 30, 942-946.

It was built as my Master's thesis at California State University, Chico.

## Overview

GraphLearnR learns the structure of a graph (its Laplacian matrix) from noisy signal observations, while accounting for external regressors that also drive the signals. This is useful for:

- **Graph signal processing** with missing data
- **Network inference** from observational data
- **Spatial data analysis** with covariates
- **Time series on graphs** with external variables

### Features

- Learn graph structure from signals with external regressors
- Handle missing data through signal estimation
- Evaluate learned graphs with precision, recall, F-score and normalized mutual information (NMI)
- Generate synthetic graphs and signals (Barabási-Albert, Erdős-Rényi) for benchmarking
- Tested on real sensor networks (NOAA weather stations, soil moisture sensors)

## Installation

```r
# install.packages("devtools")
devtools::install_github("Shivam9927/GraphLearnR")
```

Or from source:

```bash
git clone https://github.com/Shivam9927/GraphLearnR.git
```

```r
devtools::install("GraphLearnR")
```

## Quick start

```r
library(GraphLearnR)
library(igraph)

set.seed(42)

# Generate a synthetic graph
G <- sample_pa(20, m = 1, directed = FALSE)

# Generate signals on the graph with external regressors
sim <- generate_signals(G, num_signals = 100, beta = c(1, 1), seed = 42)

# Learn the graph structure
fit <- glreg(
  X_noisy = sim$S_noisy,
  P = sim$P,
  w1 = 0.001,
  w2 = 0.99
)

# Compare the learned graph with the true graph
metrics <- evaluate_graph(sim$L, fit$L)
metrics
```

`evaluate_graph()` prints precision, recall, F-score and NMI for the learned edges. Results depend on the graph, the seed and the hyperparameters `w1` and `w2`.

## Core functions

### `glreg()`

Main function: learns the graph Laplacian with external regressors.

```r
fit <- glreg(
  X_noisy,              # Observed signals (N x num_signals matrix)
  P,                    # List of regressor matrices
  w1,                   # Smoothness penalty
  w2,                   # Sparsity penalty
  max_iter = 100,       # Maximum iterations
  threshold = 0.001,    # Edge pruning threshold
  lambda_reg = 1e-6,    # Numerical regularization
  tol = 1e-4,           # Convergence tolerance
  verbose = FALSE       # Print progress
)
```

Returns an object of class `glreg` with:

- `L`: learned (pruned) Laplacian matrix
- `L_unpruned`: Laplacian before edge pruning
- `hatS`: estimated signals
- `beta`: regression coefficients
- `objective`: objective value at each iteration

`print()` and `summary()` methods are included.

### `generate_signals()`

Generates synthetic signals on a graph with external regressors.

```r
sim <- generate_signals(
  graph,                # igraph object
  num_signals,          # Number of signal observations
  beta,                 # Regression coefficients
  P = NULL,             # Regressor matrices (generated if NULL)
  mu = 0,               # Noise mean
  sigma = 0.5,          # Noise standard deviation
  regressor_mean = 0,   # Mean of generated regressors
  regressor_sd = 10,    # SD of generated regressors
  seed = NULL           # Random seed
)
```

Returns `S_clean`, `S_noisy`, `P`, the true normalized Laplacian `L`, and `beta`. `generate_signals_with_distribution()` does the same with other noise distributions.

### `evaluate_graph()`

Compares a learned Laplacian with the true one.

```r
metrics <- evaluate_graph(L_true, L_learned)
```

Returns precision, recall, F-score, NMI, and the counts of true and false positive and negative edges. `compare_methods()` compares several learned graphs at once, and `evaluate_glreg()` evaluates a `glreg` fit directly.

## Real-world application

The package was applied to temperature data from 20 NOAA weather stations, using daily maximum temperature as the signal and elevation and minimum temperature as regressors. The analysis scripts, data and figures are in [`Weather_data/`](Weather_data/).

In the thesis, GLReg improved network recovery by about 33 to 35% compared with graph learning without regressors.

## Algorithm

GLReg solves:

```
minimize  ||hatS - S||^2 + w1 * tr(hatS^T L hatS) + w2 * ||L||_F^2
subject to:
  tr(L) = N
  L_ij <= 0 for i != j
  L_ii >= 0
  L * 1 = 0
  hatS = P * beta + smooth graph signal
```

where `S` is the observed signal, `hatS` the estimated signal, `L` the graph Laplacian, `P` the external regressors, `beta` the regression coefficients, and `w1`, `w2` the hyperparameters.

The algorithm alternates between three steps until convergence:

1. Update `L` (convex optimization with CVXR)
2. Update `hatS` (closed-form solution)
3. Update `beta` (closed-form solution)

## Requirements

- R 3.5.0 or later
- CVXR 1.8.0 or later (convex optimization)
- igraph 1.2.0 or later (graph operations)
- Matrix

Optional: `mice` (missing data imputation), `ggplot2` (plots).

## Citation

```bibtex
@mastersthesis{pawar2026graphlearnr,
  title  = {GraphLearnR: An R Package for Powerful and Accessible Graph Learning},
  author = {Pawar, Shivam},
  year   = {2026},
  school = {California State University, Chico},
  type   = {Master's Thesis}
}

@article{guo2023graph,
  title   = {Graph Learning From Signals With Smoothness Superimposed by Regressors},
  author  = {Guo, Jane and Moses, Samuel and Wang, Zhihua},
  journal = {IEEE Signal Processing Letters},
  volume  = {30},
  pages   = {942--946},
  year    = {2023}
}
```

## Acknowledgments

- **Dr. Jane Guo** (California State University, Chico): thesis advisor and GLReg algorithm developer
- **Christopher Duran**: Python implementation reference
- **Prof. Samuel Moses**: methodology development

## Author

**Shivam Pawar**  
MS Data Science and Analytics, California State University, Chico  
GitHub: [@Shivam9927](https://github.com/Shivam9927)

Version 0.1.0 (May 2026). MIT License, see [LICENSE](LICENSE).
