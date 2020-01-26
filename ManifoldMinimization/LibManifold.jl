using Distributed
using JuMP, Ipopt ## Notice JuMP v0.18 is used
using Dates

function CentralizedSolver_v018(np::Int64)
    println("Centralized Ipopt with np=$np")
    m = Model(solver=IpoptSolver(linear_solver="MA27"))
    ## initial values
    startx = [0.2 for i in 1:np]
    starty = [0.3 for i in 1:np]
    startz = [0.1 for i in 1:np]
    if np == 300 ## Ipopt encounters restoration failure for default inital values
        startx = [0.5 for i in 1:np]
        starty = [0.5 for i in 1:np]
        startz = [0.4 for i in 1:np]
    end
    ## add variables
    @variable(m, x[i in 1:np], lowerbound = -1.0, upperbound = 1.0, start = startx[i])
    @variable(m, y[i in 1:np], lowerbound = -1.0, upperbound = 1.0, start = starty[i])
    @variable(m, z[i in 1:np], lowerbound = -1.0, upperbound = 1.0, start = startz[i])

    # add constraints
    @NLconstraint(m, [i=1:np], x[i]^2 +y[i]^2+z[i]^2 == 1)
    # add objevtive
    @NLobjective(m, Min,
    sum(sum(1/sqrt(2-2*x[i]*x[j]-2*y[i]*y[j]-2*z[i]*z[j]) for j in i+1:np) for i in 1:np-1))
    return m
end

## helper functions
function GenerateTotalIdx(k::Int64, np::Int64)
    Dict_total_idx = Dict()
    for k in 1:3
        if k == 1
            self_b = 1
            self_e = trunc(Int, np/3)
            other_b = trunc(Int, np/3+1)
            other_e = trunc(Int, 2*np/3)
        elseif k == 2
            self_b = trunc(Int, np/3+1)
            self_e = trunc(Int, 2*np/3)
            other_b = trunc(Int, 2*np/3+1)
            other_e = trunc(Int, np)
        else
            self_b = trunc(Int, 2*np/3+1)
            self_e = trunc(Int, np)
            other_b = trunc(Int, 1)
            other_e = trunc(Int, np/3)
        end
        Dict_total_idx[k] = union(collect(self_b:self_e), collect(other_b:other_e))
    end
    return Dict_total_idx
end

function initialization(np::Int64, Dict_total_idx::Dict)
    lmd = []
    y = []
    z = []
    blockSol = []
    for k in 1:3
        lmdDict = Dict(i=>[0.0;0.0;0.0] for i in Dict_total_idx[k])
        yDict = Dict(i=>[0.0;0.0;0.0] for i in Dict_total_idx[k])
        zDict = Dict(i=>[0.0;0.0;0.0] for i in Dict_total_idx[k])
        solDict = Dict(i=>[0.2;0.2;0.1] for i in Dict_total_idx[k])
        push!(lmd, lmdDict)
        push!(y, yDict)
        push!(z, zDict)
        push!(blockSol, solDict)
    end
    global_copy = Dict(i=>[0.5;0.5;0.5] for i in 1:np)
    location = Dict()
    for i in 1:np
        if i <= np/3
            location[i] = [1,3]
        elseif i>2*np/3
            location[i] = [2,3]
        else
            location[i] = [1,2]
        end
    end
    return (lmd, y, z, global_copy, blockSol, location)
end

@everywhere function SubproblemConstructor_v018(k::Int64, np::Int64)
    if k == 1
        self_b = 1
        self_e = trunc(Int, np/3)
        other_b = trunc(Int, np/3+1)
        other_e = trunc(Int, 2*np/3)
    elseif k == 2
        self_b = trunc(Int, np/3+1)
        self_e = trunc(Int, 2*np/3)
        other_b = trunc(Int, 2*np/3+1)
        other_e = trunc(Int, np)
    else
        self_b = trunc(Int, 2*np/3+1)
        self_e = trunc(Int, np)
        other_b = trunc(Int, 1)
        other_e = trunc(Int, np/3)
    end
    local_idx_set = union(collect(self_b:self_e), collect(other_b:other_e))
    # construct model
    m = Model(solver=IpoptSolver(print_level=0, linear_solver="MA27"))
    @variable(m, x[local_idx_set], lowerbound = -1.0, upperbound = 1.0, start = 0.2)
    @variable(m, p[local_idx_set], lowerbound = -1.0, upperbound = 1.0, start = 0.3)
    @variable(m, q[local_idx_set], lowerbound = -1.0, upperbound = 1.0, start = 0.1)
    @variable(m, cost)
    @NLconstraint(m, [i in local_idx_set], x[i]^2 +p[i]^2+q[i]^2 == 1)

    # add objective
    couple_dict = Dict()
    for i in self_b:self_e
        if k <= 2
            temp_idx_set = []
        else
            temp_idx_set = collect(other_b:other_e)
        end
        for j in local_idx_set
            if i < j
                push!(temp_idx_set, j)
            end
        end
        couple_dict[i]=temp_idx_set
    end
    @NLconstraint(m, sum(sum(1/sqrt(2 - 2*x[i]*x[j]-2*p[i]*p[j]-2*q[i]*q[j])
        for j in couple_dict[i]) for i in self_b:self_e) <= cost)
    @objective(m, Min, cost)
    return m
end

@everywhere function PreSolve_v018()
    global m
    solve(m)
end

@everywhere function SolveSubproblem_v018(k, total_idx_set, global_copy, z_k, y_k, rho)
    global m
    expr_alm_obj = @expression(m, m[:cost])
    for i in total_idx_set
        temp = [m[:x][i]; m[:p][i]; m[:q][i]] - global_copy[i] + z_k[i]
        expr_alm_obj += y_k[i]'* temp + 0.5*rho*(temp'*temp)
    end
    @objective(m, Min, expr_alm_obj)
    solve(m)
    obj = getvalue(m[:cost])
    solDict = Dict(i=>[getvalue(m[:x][i]);
                       getvalue(m[:p][i]);
                       getvalue(m[:q][i])] for i in total_idx_set)
    return (obj, solDict)
end
################################################################################
## Blow are the corresponding functions using JuMP v0.18 where MOI is used to connect solvers.
## We observe the solution query becomes substantially slower
################################################################################

## Centralized Ipopt
# function CentralizedSolver(np::Int64)
#     println("Centralized Ipopt with np=$np")
#     m = Model(with_optimizer(Ipopt.Optimizer, linear_solver = "ma57"))
#     ## initial values
#     startx = [0.5 for i in 1:np]
#     starty = [0.5 for i in 1:np]
#     startz = [0.5 for i in 1:np]
#     ## add variables
#     @variable(m, x[i in 1:np], lower_bound = -1.0, upper_bound = 1.0, start = startx[i])
#     @variable(m, y[i in 1:np], lower_bound = -1.0, upper_bound = 1.0, start = starty[i])
#     @variable(m, z[i in 1:np], lower_bound = -1.0, upper_bound = 1.0, start = startz[i])
#     # add constraints
#     @NLconstraint(m, [i=1:np], x[i]^2 +y[i]^2+z[i]^2 == 1)
#     # add objevtive
#     @NLobjective(m, Min,
#     sum(sum(1/sqrt(2-2*x[i]*x[j]-2*y[i]*y[j]-2*z[i]*z[j]) for j in i+1:np) for i in 1:np-1))
#     return m
# end
#
## subproblem constructor
# @everywhere function SubproblemConstructor(k::Int64, np::Int64)
#     if k == 1
#         self_b = 1
#         self_e = trunc(Int, np/3)
#         other_b = trunc(Int, np/3+1)
#         other_e = trunc(Int, 2*np/3)
#     elseif k == 2
#         self_b = trunc(Int, np/3+1)
#         self_e = trunc(Int, 2*np/3)
#         other_b = trunc(Int, 2*np/3+1)
#         other_e = trunc(Int, np)
#     else
#         self_b = trunc(Int, 2*np/3+1)
#         self_e = trunc(Int, np)
#         other_b = trunc(Int, 1)
#         other_e = trunc(Int, np/3)
#     end
#     local_idx_set = union(collect(self_b:self_e), collect(other_b:other_e))
#     # construct model
#     m = Model(with_optimizer(Ipopt.Optimizer, print_level = 0))
#     @variable(m, x[local_idx_set], lower_bound = -1.0, upper_bound = 1.0, start = 0.2)
#     @variable(m, p[local_idx_set], lower_bound = -1.0, upper_bound = 1.0, start = 0.3)
#     @variable(m, q[local_idx_set], lower_bound = -1.0, upper_bound = 1.0, start = 0.1)
#     @variable(m, cost)
#     @constraint(m, [i in local_idx_set], x[i]^2 +p[i]^2+q[i]^2 == 1)
#
#     # add objective
#     couple_dict = Dict()
#     for i in self_b:self_e
#         if k <= 2
#             temp_idx_set = []
#         else
#             temp_idx_set = collect(other_b:other_e)
#         end
#         for j in local_idx_set
#             if i < j
#                 push!(temp_idx_set, j)
#             end
#         end
#         couple_dict[i]=temp_idx_set
#     end
#     @NLconstraint(m, sum(sum(1/sqrt(2 - 2*x[i]*x[j]-2*p[i]*p[j]-2*q[i]*q[j])
#         for j in couple_dict[i]) for i in self_b:self_e) <= cost)
#     @objective(m, Min, cost)
#     return m
# end
#
## Presolve with default objective
# @everywhere function PreSolve()
#     global m
#     optimize!(m)
# end
## subproblem solver
# @everywhere function SolveSubproblem(k, total_idx_set, global_copy, z_k, y_k, rho)
#     global m
#     expr_alm_obj = @expression(m, m[:cost])
#     for i in total_idx_set
#         temp = [m[:x][i]; m[:p][i]; m[:q][i]] - global_copy[i] + z_k[i]
#         expr_alm_obj += y_k[i]'* temp + 0.5*rho*(temp'*temp)
#     end
#     @objective(m, Min, expr_alm_obj)
#     JuMP.optimize!(m)
#     obj = JuMP.value(m[:cost])
#     solDict = Dict(i=>[JuMP.value(m[:x][i]);
#                        JuMP.value(m[:p][i]);
#                        JuMP.value(m[:q][i])] for i in total_idx_set)
#     return (obj, solDict)
# end
# ## JuMP model passed as an argument
# @everywhere function SolveSubproblem2(m, k, total_idx_set, global_copy, z_k, y_k, rho)
#     expr_alm_obj = @expression(m, m[:cost])
#     for i in total_idx_set
#         temp = [m[:x][i]; m[:p][i]; m[:q][i]] - global_copy[i] + z_k[i]
#         expr_alm_obj += y_k[i]'* temp + 0.5*rho*(temp'*temp)
#     end
#     @objective(m, Min, expr_alm_obj)
#     JuMP.optimize!(m)
#     obj = JuMP.value(m[:cost])
#     solDict = Dict(i=>[JuMP.value(m[:x][i]);
#                        JuMP.value(m[:p][i]);
#                        JuMP.value(m[:q][i])] for i in total_idx_set)
#     return (obj, solDict)
# end
