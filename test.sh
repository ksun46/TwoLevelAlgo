
timestamp=`date +%m%d-%H%M`

## julia
julia_cmd="alias julia=/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia"
eval $julia_cmd

## warmup
cmd0="julia -e 'include(\"warmup.jl\")'"
eval $cmd0&>warmup_${timestamp}.log

## network
cmd1="julia -e 'include(\"NetworkProblem/Test_Network.jl\")'"
eval $cmd1&>Network_${timestamp}.log

## manifold
cmd2="julia -e 'include(\"ManifoldMinimization/ParallelManifold.jl\")'"
eval $cmd2&>Manifold_${timestamp}.log

## robust pca
cmd3="julia -e 'include(\"RobustTensorPCA/RobustTensorPCA.jl\")'"
eval $cmd3&>PCA_${timestamp}.log
