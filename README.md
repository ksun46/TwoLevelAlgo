# TwoLevelAlgo
Numerical Experiments for the manuscript <em>A Two-level Distributed Algorithm for Nonconvex Constrained Optimization</em>

The Github repository is [here](https://github.com/ksun46/TwoLevelAlgo)
# Description

## Hardware
The results from `Result_KS` are generated on a 2018 MacBook Pro
* Processor: 2.6 GHz Intel Core i7
* Number of cores: 6
* Memory: 16GB

## Software
* Julia v1.1.0
* JuMP v0.18: This earlier version of JuMP is used in all related experiments, where the underlying solver interface is provided by `MathProgBase`. We notice that later versions of JuMP using `MathOptInterface` will slow down our algorithms, especiall when parallization is used.
* Ipopt v3.12.8: Linear solver MA27 is used in all related experiments. If not available on the testing environment, Ipopt can be installed from [`Ipopt.jl`](https://github.com/JuliaOpt/Ipopt.jl). However, notice that this version of Ipopt uses linear solver MUMPS, which may potentially compromise the performance on large problems. See more instructions below.

All Julia pacakges can be installed by running the script `warmup.jl` (run `Pkg.add("Ipopt")` if Ipopt is not available).

## Reproduce Results
The results obtained by the aurthor can be found in the folder `Result_KS`. To reproduce results, clone this repository, and make sure a version of Julia is properly installed on the testing environment (Julia 1.1 is used by the author and other versions of Julia haven't been tested).

### For Technical Editor (NEW)
Since Julia and all dependencies are already set up on the testing environment,
the technical editor should not worry about the Ipopt issues mentioned above.
The easiest way to reproduce all results is to simply execute the bash script file `sh test.sh`.

Alternatively, the technical editor can change current directory to the folder of this project, activate a Julia session, and execute the following lines of codes.

### Install Required Packages
Execute
```julia
include("warmup.jl")
```
The software packages specified in `Project.toml` and `Manifest.toml` will be installed. If Julia reports `julia Failed to precompile ExcelReaders`
while executing `using ExcelReaders`, try
```julia
Pkg.build("ExcelReaders")
```
and then re-run the script. Though not observed by the authors, if any other pacakage reports compilation errors, please also try to run `Pkg.build("$PACKAGE_NAME")` to see if the problem can be fixed. If some problems still occur during the precomilation phase, please locate the problematic package, try `] rm $PACKAGE_NAME` to remove the package from the project, and then readd it through `] add $PACKGAE_NAME@$VERSION`(make sure the same version is added back) and build.

A customized version of Ipopt with HSL linear solvers is assumed to be already installed, so that executing `using Ipopt` should be successful without any errors. Otherwise, there are two ways to obtain Ipopt: 1. execute `Pkg.add("Ipopt")` and then `using Ipopt`. This approach downloads Ipopt binaries from [`Ipopt.jl`](https://github.com/JuliaOpt/Ipopt.jl), and Ipopt is immediately ready to be used. However, this version of Ipopt uses linear solver MUMPS, which may not be able to reproduce results reported in the folder `Result_KS`. 2. Obtain Ipopt and HSL thrid party library according to this [documentation](https://coin-or.github.io/Ipopt/INSTALL.html). Then follow the [Custom Installation section](https://github.com/JuliaOpt/Ipopt.jl) to use customized Ipopt in Julia; alternatively, if following the custom installation does not work, one can try the following: go to, e.g., `~/.julia/packages/Ipopt/Iu7vT/deps` (`Iu7vT` is a folder that contains a version of Ipopt that will be used in this project; make sure from Julia console output that this is the version being used when running this project), open `deps.jl`, and point the variables `libipopt` and `amplexe` to the dylib and binary files of the customized Ipopt. When running `warmup.jl`, please also check the screen output to make sure Ipopt and linear solver MA27 are installed correctly.

### Network Problem
In this part, one master and up to four workers threads will be used. So a total of five threads should be available. Execute
```julia
include("NetworkProblem/Test_Network.jl")
```
This will generate `NetworkResult$(mmdd-HHMM).csv`. This script takes around 1 hr.
### Parallel Minimization over Compact Manifold
In this part, one master and three workers threads will be used. So a total of four threads should be available. Execute
```julia
include("ManifoldMinimization/ParallelManifold.jl")
```
This will generate `ParallelManifoldResult$(mmdd-HHMM).csv`. This script takes around 2 hrs.
### Robust Tensor PCA
This part is implemented in single thread. Execute
```julia
include("RobustTensorPCA/RobustTensorPCA.jl")
```
This will generate `combined_pca_err.pdf` and `combined_pca_res.pdf`. This script takes around 1 hrs.
