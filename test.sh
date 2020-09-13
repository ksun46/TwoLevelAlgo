
timestamp=`date +%m%d-%H%M`
mkdir Log_VM_${timestamp}
## julia
# julia_cmd="alias julia=/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia"
# for MPC VM
julia_cmd="alias julia=/opt/julia-1.1.0/bin/julia"
eval $julia_cmd
sleep 10
## warmup
cmd0="julia -e 'include(\"warmup.jl\")'"
eval $cmd0&>Log_VM_${timestamp}/warmup_${timestamp}.log
sleep 10
## network
cmd1="julia -e 'include(\"NetworkProblem/Test_Network.jl\")'"
eval $cmd1&>Log_VM_${timestamp}/Network_${timestamp}.log
sleep 10
## manifold
cmd2="julia -e 'include(\"ManifoldMinimization/ParallelManifold.jl\")'"
eval $cmd2&>Log_VM_${timestamp}/Manifold_${timestamp}.log
sleep 10
## robust pca
cmd3="julia -e 'include(\"RobustTensorPCA/RobustTensorPCA.jl\")'"
eval $cmd3&>Log_VM_${timestamp}/PCA_${timestamp}.log
sleep 10
