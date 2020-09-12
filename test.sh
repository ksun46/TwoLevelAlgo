julia_cmd="alias julia=/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia"
eval $julia_cmd

## warmup
cmd0="julia -e 'include(\"warmup.jl\")'"
eval $cmd0

## network
cmd1="julia -e 'include(\"NetorkProblem/Test_Network.jl\")'"
eval $cmd1

## manifold
cmd2="julia -e 'include(\"ManifoldMinimization/ParallelManifold.jl\")'"
eval $cmd2

## robust pca
cmd3="julia -e 'include(\"RobustTensorPCA/RobustTensorPCA.jl\")'"
eval $cmd3
