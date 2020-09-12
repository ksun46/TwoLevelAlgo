## This is the warm-up file
println("Start project warm up")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.status()

using Distributed
using JuMP      ## An earlier version JuMP v0.18.6 should be used for this project
using Plots
using TensorToolbox
using Dates
using MathOptInterface
using MathProgBase
using DataFrames
using ExcelReaders
using CSV

## if a customized version of Ipopt with MA27 is already installed, then directly "using Ipopt";
## otherwise run '''Pkg.add("Ipopt")''', where default linear solver MUMPS will be used,
## and perforamnce may be comproised
using Ipopt
ENV["JULIA_IPOPT_LIBRARY_PATH"] = "/Users/kaizhaosun/coinOR/Ipopt/build/lib"
ENV["JULIA_IPOPT_EXECUTABLE_PATH"] =  "/Users/kaizhaosun/coinOR/Ipopt/build/bin"
Pkg.build("Ipopt")
## test IPOPT
println("Test IPOPT in manager process...")

## JuMP 0.20
# m_test = Model(with_optimizer(Ipopt.Optimizer, linear_solver = "ma27"))
# @variable(m_test, x[i in 1:10], lower_bound = 0.0, upper_bound = 1.0)
# @constraint(m_test, sum(x) == 1.0)
# @objective(m_test, Min, sum(i * x[i] for i in 1:10))
# JuMP.optimize!(m_test)
# @assert(JuMP.termination_status(m_test)== MOI.LOCALLY_SOLVED, "Ipopt test failed. Double check solver set-up.")
# println("Ipopt works fine.")

## JuMP v0.18
m_test = Model(solver=IpoptSolver(linear_solver="MA27"))
@variable(m_test, x[i in 1:10], lowerbound = 0.0, upperbound = 1.0)
@constraint(m_test, sum(x) == 1.0)
@objective(m_test, Min, sum(i * x[i] for i in 1:10))
solve(m_test)
println("Please check Ipopt output.")
