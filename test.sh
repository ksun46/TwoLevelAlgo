julia_cmd="alias julia="
eval $julia_cmd

## network
cmd1="julia -e 'include(\"NetorkProblem/Test_Network.jl\")'"
eval $cmd1

## manifold
cmd2="julia -e 'include(\"ManifoldMinimization/ParallelManifold.jl\")'"
eval $cmd2

## robust pca
cmd3="julia -e 'include(\"RobustTensorPCA/RobustTensorPCA.jl\")'"
eval $cmd3
