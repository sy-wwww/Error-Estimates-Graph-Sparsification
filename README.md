# Error-Estimates-for-Graph-Sparsification

This repository contains MATLAB implementations of algorithms for estimating errors in graph sparsification, as described in the paper **"Empirical Error Estimates for Graph Sparsification"**. 

## Files

### 1. `bootstrap_double.m`
This file implements Algorithm 1 from the paper. It performs a double bootstrapping process to estimate errors in graph sparsification.

#### Inputs:
- `n`: Number of nodes in the graph.
- `B`: Bootstrapping sample size for the first loop.
- `C`: Bootstrapping sample size for the second loop.
- `num_func`: Number of functions to be estimated.
- `cedge`: An $N \times 4$ matrix where:
  - The first two columns represent the nodes of the sampled edges.
  - The third column contains the weights of the edges in the sparsified graph.
  - The fourth column indicates the number of times each edge is sampled.
- `y`: Parameter for graph-structured regression (optional, default = 0).
- `lambda`: Parameter for graph-structured regression (optional, default = 0)


#### Outputs:
- `er18_sd`: Error estimates for the double bootstrapping process.

---

### 2. `graph_cut.m`
This file implements Algorithm 2 from the paper. It focuses on graph cut operations to estimate errors in graph sparsification.

#### Inputs:
- `n`: Number of nodes in the graph.
- `B`: Bootstrapping sample size.
- `cedge`: An $N \times 4$ matrix where:
  - The first two columns represent the nodes of the sampled edges.
  - The third column contains the weights of the edges in the sparsified graph.
  - The fourth column indicates the number of times each edge is sampled.
- `CCS`: A sparse matrix where each row is an `n`-dimensional binary vector representing a cut in the graph.

#### Outputs:
- `er1B`: Error estimates for the bootstrapping process.


### 3. `spectral_clustering.m`
This file implements Algorithm 2 for spectral clustering from the paper. It performs spectral clustering on the graph.

#### Inputs:
- `n`: Number of nodes in the graph.
- `B`: Bootstrapping sample size for the first loop.
- `cedge`: An $N \times 4$ matrix where:
  - The first two columns represent the nodes of the sampled edges.
  - The third column contains the weights of the edges in the sparsified graph.
  - The fourth column indicates the number of times each edge is sampled.
- `r`: A number that the user believes is safely above the correct number of clusters.

#### Outputs:
- `er1B`:  Error estimates for spectral clustering, i.e. $\xi_b^*= \max_{2 \leq j \leq r} |\lambda_j(\hat{L}^*)/\lambda_j(\hat L)-1|$.
