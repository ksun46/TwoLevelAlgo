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
* JuMP v0.18: This earlier version of JuMP is used in all related experiments, where the underlying solver interface is provided by `MathProgBase`. We notice that later versions of JuMP using `MathOptInterface` will slow down our algorithms, especiall when parallization is used.
* Ipopt v3.12.8: Linear solver MA27 is used in all related experiments. If not available on the testing environment, Ipopt can be installed from [`Ipopt.jl`](https://github.com/JuliaOpt/Ipopt.jl), where custom installation is also discussed. See also [here](https://coin-or.github.io/Ipopt/INSTALL.html) on how to obtain the HSL thrid party library.
Notice that the Ipopt uses default linear solver MUMPS, which may potentially compromise the performance on large problems.

All other Julia pacakges can be installed by running the script `warmup.jl`

## Reproduce Results 
Change current directory to the folder of this project, and activate a Julia session.
### Install Required Packages
```
julia>include("warmup.jl")
```
The software packages specified in `Project.toml` and `Manifest.toml` will be installed. Please also check the screen output to make sure Ipopt is installed correctly.
### Network Problem
In this part, one master and up to four workers threads will be used. So a total of five threads should be available.
```
julia>include("NetorkProblem/Test_Network.jl")
```
This will generate `NetworkResult$(mmdd-HHMM).csv`
### Parallel Minimization over Compact Manifold
In this part, one master and three workers threads will be used. So a total of four threads should be available.
```
julia>include("ManifoldMinimization/ParallelManifold.jl")
```
This will generate `ParallelManifoldResult$(mmdd-HHMM).csv`
### Robust Tensor PCA
This part is implemented in single thread.
```
julia>include("RobustTensorPCA/Test_Network.jl")
```
This will generate `combined_pca_err.pdf` and `combined_pca_res.pdf`
