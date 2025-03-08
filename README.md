# Error-Estimates-for-Graph-Sparsification

This repository contains MATLAB implementations of algorithms for estimating errors in graph sparsification, as described in the paper **`Empirical Error Estimates for Graph Sparsification'**. 

## Files

### 1. `bootstrap_double.m`
This file implements **Algorithm 1** from the paper. It performs a double bootstrapping process to estimate errors in graph sparsification.

#### Inputs:
- `n`: Number of nodes in the graph.
- `B`: Bootstrapping sample size for the first loop.
- `C`: Bootstrapping sample size for the second loop.
- `num_func`: Number of functions to be estimated.
- `cedge`: An `N×4` matrix where:
  - The first two columns represent the nodes of the sampled edges.
  - The third column contains the weights of the edges in the sparsified graph.
  - The fourth column indicates the number of times each edge is sampled.
- `y`: Parameter for graph-structured regression (optional, default = 0).
- `lambda`: Parameter for graph-structured regression (optional, default = 0)


#### Outputs:
- `er18_sd`: Error estimates for the double bootstrapping process.

#### Key Steps:
1. Constructs a sparse matrix `S` from the input edge data.
2. Performs bootstrapping to estimate errors and their standard deviations.

---

### 2. `graph_cut.m`
This file implements **Algorithm 2** from the paper. It focuses on graph cut operations to estimate errors in graph sparsification.

#### Inputs:
- `n`: Number of nodes in the graph.
- `B`: Bootstrapping sample size.
- `cedge`: An `N×4` matrix where:
  - The first two columns represent the nodes of the sampled edges.
  - The third column contains the weights of the edges in the sparsified graph.
  - The fourth column indicates the number of times each edge is sampled.
- `CCS`: A sparse matrix where each row is an `n`-dimensional binary vector representing a cut in the graph.

#### Outputs:
- `er1B`: Error estimates for the bootstrapping process.
