using Distributed
## Creare additional worker processes
time_compilation_start = time()
println("Start compilation...")
Max_threads = Sys.CPU_THREADS
println("  Max number of threads is $Max_threads, 4 of which will be used.")
available_threads = 4 #trunc(Int, Sys.CPU_THREADS*0.75)
if Distributed.nprocs() < available_threads
    println("  Add worker processes ...")
    Distributed.addprocs(available_threads - Distributed.nprocs())
end
println("  Manager id is: $(myid())")
workersID = Distributed.workers()
println("  Workers: $(workersID)")
## Create map from zone id to worker id
MapZoneToWorker = Dict(k=>workersID[k] for k in 1:available_threads-1)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using JuMP, Ipopt
using Dates
# @everywhere using MathOptInterface
include("LibManifold.jl")
import LinearAlgebra:norm
println("Compilation finished ($(time() - time_compilation_start)).")


function ParallelManifold(np::Int64, bool_dual_update::Bool)
    ## check if num of workers is enough
    if 3 > available_threads - 1
        println("Warning: This test case requires 4 workers.")
        @assert(4<= available_threads, "Need more workers!")
    end
    if bool_dual_update
        println("Outer-level ALM for np=$np")
    else
        println("Outer-level penalty method for np=$np")
    end
    println("Start algorithm initialization...")
    time_init_start = time()

    ## create and send model to worker thread
    println("  Send subproblem models to worker processes ...")
    global m
    for k in 1:3
        println("    Construct subproblem model $k in the manager process...")
        m = SubproblemConstructor_v018(k, np)
        ## for JuMP 0.20
        # m = SubproblemConstructor(k, np)
        worker_id = MapZoneToWorker[k]
        println("    Send zone $k subproblem model to worker $worker_id...")
        remotecall_fetch(()->m, worker_id)
    end
    ## Presolve to warmup IPOPT in each worker process
    println("  Presolve to warm-up IPOPT in each worker process ...")
    @sync for k in 1:3
        @async begin
            worker_id = MapZoneToWorker[k]
            remotecall_fetch(PreSolve_v018, worker_id)
            # for JuMP 0.20
            # remotecall_fetch(PreSolve, worker_id)
        end
    end

    ## create necessary data structures
    Dict_total_idx = GenerateTotalIdx(3, np)
    (lmd, y, z, global_copy, blockSol, location) = initialization(np, Dict_total_idx)
    println("Algorithm initialization finished ($(time() - time_init_start)).")

    ## some algorithmic parameters
    maxIter = 500
    theta = 0.5 ## residual decrease factor
    gamma = 2.0 ## penalty increase factor
    epi =1.0e-6 ## tolerance
    if np < 100 ## inital penalty beta
        beta = 100.0
    elseif np < 200
        beta = 200.0
    else
        beta = 500.0
    end
    rho = 2.0 * beta     ## inner ADMM penalty
    lmd_bd = 1.0e6       ## bound on outer dual
    iterCount = 1
    almCount = 1
    slack_z_prev = 0.0
    obj_cost = 0.0
    res = 0.0
    time_loop_start = time()
    while iterCount <= maxIter
        println("Iteration: ", iterCount)
        objDict = [0.0;0.0;0.0]
        ## update the first block in parallel
        @sync for k in 1:3
            @async begin
                (objDict[k], blockSol[k])=remotecall_fetch(SolveSubproblem_v018, MapZoneToWorker[k],
                    k, Dict_total_idx[k], global_copy, z[k], y[k], rho)
                ## For JuMP 0.20
                # (objDict[k], blockSol[k])=remotecall_fetch(SolveSubproblem, MapZoneToWorker[k],
                #     k, Dict_total_idx[k], global_copy, z[k], y[k], rho)
            end
        end
        obj_cost = sum(objDict)
        ## update the second block
        for i in 1:np
        global_copy[i] = (rho*(sum(blockSol[k][i]+z[k][i] for k in location[i])) +
            sum(y[k][i] for k in location[i]))/(length(location[i])*rho)
        end

        ## update the third block and dual variable
        re2_list = []
        re3_list = []
        z_list = []
        for k in 1:3
            for i in Dict_total_idx[k]
                z[k][i] =(rho*(global_copy[i] - blockSol[k][i])-lmd[k][i]-y[k][i])/(beta+rho)
                re2_temp = blockSol[k][i]-global_copy[i]
                re3_temp = blockSol[k][i]-global_copy[i]+ z[k][i]
                y[k][i] = y[k][i] + rho * re3_temp
                push!(re2_list, re2_temp)
                push!(re3_list, re3_temp)
                push!(z_list, z[k][i])
            end
        end
        println("   Current ||Ax+Bx_bar+z||: $(norm(re2_list))")
        println("   Current   ||Ax+Bx_bar||: $(norm(re3_list))")
        println("   Current           ||z||: $(norm(z_list))")
        ## check inner termination
        if norm(re3_list) <= sqrt(3*np)/(2500*almCount)
            println("   ADMM terminates at iteration: ", iterCount)
            println("   ALM iteration $almCount finished")
            println("")
            ## check outer termination
            if norm(re2_list) <= 1e-6 * sqrt(3 * np)
                println("Converged successfully. ALM Stopping Criteria Met at Iteration $almCount")
                res = norm(re2_list)
                break
            end
            ## update outer dual
            if bool_dual_update
                for k in 1:3
                    for i in Dict_total_idx[k]
                        lmd[k][i] = max.(min.(lmd[k][i]+beta*z[k][i], lmd_bd), -lmd_bd)
                    end
                end
            end
            ## update penalty
            if norm(z_list) > theta * slack_z_prev
                beta = gamma * beta
                rho = 2.0 * beta
            end
            ## initalize next ADMM
            for k in 1:3
                for i in Dict_total_idx[k]
                    z[k][i] = [0.0; 0.0; 0.0]
                    y[k][i] = -lmd[k][i]
                end
            end
            almCount += 1
        end
        slack_z_prev = norm(z_list)
        iterCount += 1
    end
    duration = time() - time_loop_start
    println("")
    println("Number of inner ADMM iterations: $iterCount/$maxIter")
    println("Number of outer ALM  iterations: $almCount")
    println("Objevtive value at termination: $obj_cost")
    println("    Penalty_beta at termination: $beta")
    println("The Two-level Algorithm finished in $duration seconds.")
    solDict = Dict("Time"=>duration, "Obj"=>obj_cost, "Inner"=>iterCount,
                   "Outer"=>almCount, "Res"=>res)
    return solDict
end

function WriteManifoldOutput(solDict, np_list)
    now = Dates.format(Dates.now(), "mmdd-HHMM")
    output_dir = "ParallelManifoldResult$now.csv"
    open(output_dir, "w") do f
        write(f, "np, Central_Obj, Time, Method, Outer, Inner, Res, Obj,Gap, Time\n")
        for np in np_list
            central_obj = solDict[np]["central"]["obj"]
            central_time = solDict[np]["central"]["time"]

            tl_inner = solDict[np]["ALM"]["Inner"]
            tl_outer = solDict[np]["ALM"]["Outer"]
            tl_res = solDict[np]["ALM"]["Res"]
            tl_obj = solDict[np]["ALM"]["Obj"]
            tl_time = solDict[np]["ALM"]["Time"]
            tl_gap = (tl_obj - central_obj)/tl_obj
            p_inner = solDict[np]["Penalty"]["Inner"]
            p_outer = solDict[np]["Penalty"]["Outer"]
            p_res = solDict[np]["Penalty"]["Res"]
            p_obj = solDict[np]["Penalty"]["Obj"]
            p_time = solDict[np]["Penalty"]["Time"]
            p_gap = (p_obj - central_obj)/p_obj
            write(f, "$np, $central_obj, $central_time, Proposed, $tl_inner, $tl_outer, $tl_res, $tl_obj, $tl_gap, $tl_time \n", )
            write(f, " , , , Penalty, $p_inner, $p_outer, $p_res, $p_obj,$p_gap, $p_time \n", )
        end
    end
end

## run test
np_list = [60, 90, 120, 180, 240, 300]
solDict = Dict()
for np in np_list
    tempDict=Dict()
    m_central= CentralizedSolver_v018(np)
    start = time()
    solve(m_central)
    time_central = time()-start
    tempDict["central"] = Dict("time"=>time_central, "obj"=>getobjectivevalue(m_central))
    solDict_penalty = ParallelManifold(np, false)
    tempDict["Penalty"] = solDict_penalty
    solDict_alm = ParallelManifold(np, true)
    tempDict["ALM"] = solDict_alm
    solDict[np] = tempDict
end

## Ouput solution
WriteManifoldOutput(solDict, np_list)
