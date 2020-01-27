using Distributed
## Creare additional worker processes
time_compilation_start = time()
println("Start compilation...")
Max_threads = Sys.CPU_THREADS
println("  Max number of threads is $Max_threads, 4 of which will be used.")
available_threads = 5 #trunc(Int, Sys.CPU_THREADS*0.75)
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
include("LibNetwork.jl")


function PDD_Network(case, num_partition)
    data_network, Network = parser(case, num_partition)
    Dict_Nbr, List_DiLine = GetNetworkInfo(data_network)
    partition, Dict_Location_Line, Dict_Location = GeneratePartitionInfo(Network, num_partition, data_network)
    println("Start algorithm initialization...")
    time_init_start = time()
    println("  Send subproblem models to worker processes ...")
    global m
    for k in 1:num_partition
        println("    Construct subproblem model $k in the manager process...")
        m = NetworkSubproblemConstructor(partition[k], data_network, Dict_Nbr)
        worker_id = MapZoneToWorker[k]
        println("    Send zone $k subproblem model to worker $worker_id...")
        remotecall_fetch(()->m, worker_id)
    end
    ## Presolve to warmup IPOPT in each worker process
    println("  Presolve to warm-up IPOPT in each worker process ...")
    @sync for k in 1:num_partition
        @async begin
            worker_id = MapZoneToWorker[k]
            remotecall_fetch(PreSolve_Network, worker_id)
        end
    end
    Dual_lmd, Dual_y, Slack_z, BlockSol, Global_copy=Initialization(partition, Dict_Location_Line, Dict_Location)
    println("Algorithm initialization finished ($(time() - time_init_start)).")

    iterCount = 1
    almCount = 1
    epi = 1e-5
    rho = 2000.0
    theta = 0.75
    gamma = 1.5
    MaxIter = 1000
    time_loop_start = time()
    gen_cost = 0.0
    res = 0.0
    re2_prev = 0.0
    eta = 100.0
    L_prev = 1.0e12
    dim_couple = sum(length(Dict_Location[i]) for i in keys(Dict_Location)) +
        2 * length(keys(Dict_Location_Line))
    while iterCount <= MaxIter
        println("Inner BSUM iteration $iterCount starts...")
        ## Sample random block index
        block_type, block_idx = SampleAndUpdateFirstBlock_PDD!(Dict_Location,
            Dict_Location_Line, BlockSol, Global_copy, Dual_y, rho)
        ## Solve nonlinear subproblems
        obj_list = [0.0 for k in 1:num_partition]
        @sync for k in 1:num_partition
            @async begin
                BlockSol[k], obj_list[k] = remotecall_fetch(SolveSubproblem_PDD,
                    MapZoneToWorker[k], partition[k], Global_copy, Dual_y[k],rho)
            end
        end
        gen_cost = sum(obj_list)
        ## Update second block
        UpdateSecondBlock_PDD!(Global_copy, BlockSol, Dual_y, Dict_Location_Line,
            Dict_Location, rho, block_idx)
        re2, L_obj = Calculate_metrics!(BlockSol, Global_copy, Dual_y, gen_cost,rho)
        println("   Current Generation Cost: $gen_cost")
        println("   Current   ||Ax+Bx_bar||: $(re2)")
        res = re2
        ## chcek inner stop criteria
        if abs(L_obj-L_prev)/ abs(L_prev) <= max(epi,1.0e-3 * (2/3)^(almCount))
            println("   BSUM terminates at iteration: ", iterCount)
            println("   ALM iteration $almCount finished")
            println("   Current rho is $rho")
            println("")
            ## check outer stop criteria
            if re2 <= epi * sqrt(dim_couple)
                println("Converged successfully. ALM Stopping Criteria Met at Iteration $almCount")
                break
            end
            eta = theta * min(eta, re2_prev)
            if re2 <= eta
                UpdateDual_PDD!(BlockSol, Global_copy, Dual_y,rho)
            else
                rho = rho * gamma
            end
            re2_prev = re2
            almCount += 1
        end
        iterCount +=1
        L_prev = L_obj
    end
    duration = time() - time_loop_start
    iterCount = min(iterCount, MaxIter)
    println("")
    println("Number of inner BSUM iterations: $iterCount")
    println("Number of outer ALM  iterations: $almCount")
    println("Objevtive value at termination: $gen_cost")
    println("    Penalty rho at termination: $rho")
    println("The PDD Algorithm finished in $duration seconds.")
    return Dict("Inner"=>iterCount, "Outer"=>almCount, "Res"=>res,
                "Cost"=>gen_cost, "Time"=>duration)
end

#
# case="case300"
# num_partition = 3
# PDD_Network(case, num_partition)
