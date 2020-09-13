## time stamp
timestamp=`date +%m%d-%H%M`
mkdir Log_VM_${timestamp}

# for MPC VM
julia_cmd="alias julia=/opt/julia-1.1.0/bin/julia"
eval $julia_cmd

## warmup
julia warmup.jl>Log_VM_${timestamp}/warmup_${timestamp}.log

## network
julia -e 'include(\"NetworkProblem/Test_Network.jl\")'>Log_VM_${timestamp}/Network_${timestamp}.log
## manifold
julia -e 'include(\"ManifoldMinimization/ParallelManifold.jl\")'>Log_VM_${timestamp}/Manifold_${timestamp}.log
## pca
julia -e 'include(\"RobustTensorPCA/RobustTensorPCA.jl\")'>Log_VM_${timestamp}/PCA_${timestamp}.log


## julia cmd
# julia_cmd="alias julia=/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia"
# eval $julia_cmd
## network
# cmd1="julia -e 'include(\"NetworkProblem/Test_Network.jl\")'"
# eval $cmd1>Log_VM_${timestamp}/Network_${timestamp}.log
## manifold
# cmd2="julia -e 'include(\"ManifoldMinimization/ParallelManifold.jl\")'"
# eval $cmd2>Log_VM_${timestamp}/Manifold_${timestamp}.log
## robust pca
# cmd3="julia -e 'include(\"RobustTensorPCA/RobustTensorPCA.jl\")'"
# eval $cmd3>Log_VM_${timestamp}/PCA_${timestamp}.log
