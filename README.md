# TwoLevelAlgo
Numerical Experimeents for the manuscript <em>A Two-level Distributed Algorithm for Nonconvex Constrained Optimization</em>

# Description

## Hardware

## Software
* Julia v1.1.0
* JuMP v0.18: This earlier version of JuMP is used in all related experiments, where the underlying solver interface is provided by `MathProBase`. We notice 
* Ipopt 3.12.8: Linear solver MA27 is used in all related experiments. If not available on the testing environment, Ipopt can be installed from [`Ipopt.jl`](https://github.com/JuliaOpt/Ipopt.jl), where custom installation is also discussed. Notice that the Ipopt uses default linear solver MUMPS, which may potentially compromise the performance on large problems.

All other Julia pacakges can be installed by running the script `warmup.jl`

## Reproduce Results 
### Install Required Packages
```
include("warmup.jl")

```
### Network Problem
```
include("NetorkProblem/Test_Network.jl")

```

### Parallel Minimization over Compact Manifold
```
include("ManifoldMinimization/ParallelManifold")

```

### Robust Tensor PCA
```
include("RobustTensorPCA/Test_Network.jl")

```
