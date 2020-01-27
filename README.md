# TwoLevelAlgo
Numerical Experimeents for the manuscript <em>A Two-level Distributed Algorithm for Nonconvex Constrained Optimization</em>

# Description

## Hardware
The results from `Result_kS` are generated on a 2018 MacBook Pro 
* Processor: 2.6 GHz Intel Core i7
* Number of cores: 6
* Memory: 16GB

## Software
* Julia v1.1.0
* JuMP v0.18: This earlier version of JuMP is used in all related experiments, where the underlying solver interface is provided by `MathProBase`. We notice 
* Ipopt v3.12.8: Linear solver MA27 is used in all related experiments. If not available on the testing environment, Ipopt can be installed from [`Ipopt.jl`](https://github.com/JuliaOpt/Ipopt.jl), where custom installation is also discussed. Notice that the Ipopt uses default linear solver MUMPS, which may potentially compromise the performance on large problems.

All other Julia pacakges can be installed by running the script `warmup.jl`

## Reproduce Results 
Change current directory to the folder of this project, and activate a Julia session.
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
