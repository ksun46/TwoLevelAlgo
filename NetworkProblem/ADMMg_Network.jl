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


function ADMMg_Network(case, num_partition)
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
    Dict_total_gen = Dict()
    for k in 1:num_partition
        total_gen = []
        for i in partition[k]["SelfGen"]
            for g in data_network["Gen"]["GenIdx"][i]
                push!(total_gen, g)
            end
        end
        Dict_total_gen[k] = total_gen
    end
    Dual_lmd, Dual_y, Slack_z, BlockSol, Global_copy=Initialization(partition, Dict_Location_Line, Dict_Location)
    println("Algorithm initialization finished ($(time() - time_init_start)).")
    iterCount = 1
    almCount = 1
    epi= 1e-5
    beta = 1/(epi)^2 ## last block penalty
    rho = 3*beta     ## ADMM penalty
    h = 0.01/epi     ## proximal coefficient
    step_size = 1/(rho)
    MaxIter = 1000
    time_loop_start = time()
    gen_cost = 0.0
    res = 0.0
    rez_prev = 0.0
    dim_couple = sum(length(Dict_Location[i]) for i in keys(Dict_Location)) +
        2 * length(keys(Dict_Location_Line))
    while iterCount <= MaxIter
        println("ADMMg iteration $iterCount starts.")
        obj_list = [0.0 for k in 1:num_partition]
        @sync for k in 1:num_partition
            @async begin
                BlockSol[k], obj_list[k] = remotecall_fetch(SolveSubproblem_ADMMg,
                    MapZoneToWorker[k], partition[k], Global_copy, Dict_total_gen[k],
                    Dual_y[k], Slack_z[k], rho, h)
            end
        end
        gen_cost = sum(obj_list)
        ## second block update
        UpdateSecondBlock_ADMMg!(Global_copy, BlockSol, Slack_z, Dual_y, Dict_Location_Line, Dict_Location, rho, h)
        ## thrid block update
        UpdateThirdBlock_ADMMg!(Slack_z, BlockSol, Global_copy, Dual_y, rho, beta, step_size)
        ## dual update
        re2, re3, rez = UpdateDual_y!(Dual_y, BlockSol, Global_copy,Slack_z, rho)
        println("   Current Generation Cost: $gen_cost")
        println("   Current ||Ax+Bx_bar+z||: $(re3)")
        println("   Current   ||Ax+Bx_bar||: $(re2)")
        println("   Current           ||z||: $(rez)")
        res = re2
        ## chcek inner stop criteria
        if re2 <= epi * sqrt(dim_couple)
            println("Converged successfully. ADMM_g Stopping Criteria Met at Iteration $iterCount")
            break
        end
        iterCount +=1
    end
    iterCount = min(iterCount, MaxIter)
    duration = time() - time_loop_start
    println("")
    println("Number of inner ADMM iterations: $iterCount")
    println("Objevtive value at termination: $gen_cost")
    println("The ADMM-g Algorithm finished in $duration seconds.")
    return Dict("Inner"=>iterCount, "Res"=>res, "Cost"=>gen_cost, "Time"=>duration)
end

#
# case="case300"
# num_partition = 2
# ADMMg_Network(case, num_partition)
