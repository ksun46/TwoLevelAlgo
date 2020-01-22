## This is the warm-up file
println("Start project warm up")
time_start = time()
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.status()

using Distributed
using JuMP
using Plots
using TensorToolbox
using Dates
using MathOptInterface
## if a customized version of Ipopt is already installed, then directly "using Ipopt"
## otherwise run '''Pkg.add("Ipopt")''', where default linear solver MUMPS will be used
using Ipopt

## test IPOPT
println("Test IPOPT in manager process...")
m_test = Model(with_optimizer(Ipopt.Optimizer, linear_solver = "ma57"))
@variable(m_test, x[i in 1:10], lower_bound = 0.0, upper_bound = 1.0)
@constraint(m_test, sum(x) == 1.0)
@objective(m_test, Min, sum(i * x[i] for i in 1:10))
JuMP.optimize!(m_test)
@assert(JuMP.termination_status(m_test)== MOI.LOCALLY_SOLVED, "Ipopt test failed. Double check solver set-up.")
println("Ipopt works fine.")
println("Project warm-up finished ($(time() - time_start)).")
