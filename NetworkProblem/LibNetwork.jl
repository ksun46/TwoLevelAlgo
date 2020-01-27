using Distributed
using DataFrames, ExcelReaders, CSV
using Ipopt, JuMP
const mBase = 100.0
import LinearAlgebra:norm
## Parser
function GenerateLineDict(fromBus, toBus, flow_max)
    Dict_Line = Dict()
    for i in 1:length(fromBus)
        temp = (min(fromBus[i],toBus[i]), max(fromBus[i],toBus[i]))
        Dict_Line[temp] = Float64(flow_max[i])
    end
    return Dict_Line
end

function GenerateGenDict(gen_idx)
    Dict_GenIdx = Dict(i=>[] for i in unique(gen_idx))
    for i in 1:length(gen_idx)
        push!(Dict_GenIdx[gen_idx[i]], i)
    end
    return Dict_GenIdx
end

function GenerateYMatrix(f,t,G_val,B_val)
    G = Dict()
    B = Dict()
    for l in 1:length(f)
        i = f[l]
        j = t[l]
        G[(i,j)] = G_val[l]
        B[(i,j)] = B_val[l]
    end
    return (G,B)
end

function parser(case, num_partition)
    ## Read raw data
    dir1 = string("NetworkProblem/data/$case/",string(case),".xlsx")
    dir2 = string("NetworkProblem/data/$case/",string(case),"partition.xlsx")
    line = readxlsheet(dir1, "branch", skipstartrows=1)
    bus = readxlsheet(dir1, "bus", skipstartrows = 1)
    gen = readxlsheet(dir1, "gen", skipstartrows = 1)
    Y = readxlsheet(dir1,"Y", skipstartrows = 1)
    Network = readxlsheet(dir2, string("Sheet",num_partition-1),skipstartrows=1)

    ## Bus
    busNumber = [Int64(item) for item in bus[:,1]]
    load_P = bus[:,3]/mBase
    load_Q = bus[:,4]/mBase
    v_max = bus[:,9]
    v_min = bus[:,10]
    Dict_Bus = Dict("BusIdx"=>busNumber,
                    "Load_P"=>load_P, "Load_Q"=>load_Q,
                    "v_max"=>v_max,   "v_min"=>v_min)
    ## Branch
    fromBus = [Int64(item) for item in line[:,1]]
    toBus = [Int64(item) for item in line[:,2]]
    flow_max = line[:,3]/100
    Dict_Line = GenerateLineDict(fromBus, toBus, flow_max)

    ## Generator
    gen_idx = [Int64(item) for item in gen[:,1]]
    p_max = gen[:,4]/mBase
    p_min = gen[:,5]/mBase
    q_max = gen[:,2]/mBase
    q_min = gen[:,3]/mBase
    c2 = gen[:,7]
    c1 = gen[:,8]
    c0 = gen[:,9]
    Dict_GenIdx = GenerateGenDict(gen_idx)
    Dict_Gen = Dict("GenIdx"=> Dict_GenIdx, "p_max"=>p_max, "p_min"=>p_min,
                    "q_max"=>q_max, "q_min"=>q_min, "c2"=>c2, "c1"=>c1, "c0"=>c0)

    ## Y bus
    f = [Int64(item) for item in Y[:,1]]
    t = [Int64(item) for item in Y[:,2]]
    G_val = Y[:,3]
    B_val = Y[:,4]
    (G, B) = GenerateYMatrix(f,t,G_val, B_val)
    data_network = Dict("Bus"=>Dict_Bus, "Line"=>Dict_Line, "Gen"=>Dict_Gen,
                        "G"=>G, "B"=>B)
    return data_network, Network
end

## Get additional information
function GetNetworkInfo(data_network)
    Dict_Bus  = data_network["Bus"]
    Dict_Line = data_network["Line"]
    Dict_Gen  = data_network["Gen"]
    ## find neighbors of each bus and directed lines
    Dict_Nbr = Dict(i=>[] for i in 1:length(Dict_Bus["BusIdx"]))
    List_DiLine = []
    for (i,j) in keys(Dict_Line)
        push!(Dict_Nbr[i], j)
        push!(Dict_Nbr[j], i)
        push!(List_DiLine, (i,j))
        push!(List_DiLine, (j,i))
    end
    return Dict_Nbr, List_DiLine
end

function isFloat(x)
    if typeof(x) == Float64
        return true
    else
        return false
    end
end

## Generate partition Info
function GeneratePartitionInfo(Network, num_partition, data_network)
    Dict_Bus  = data_network["Bus"]
    Dict_Line = data_network["Line"]
    Dict_Gen  = data_network["Gen"]
    partition = Dict(k =>Dict() for k in 1:num_partition)
    Dict_BusToZone = Dict()
    for k in 1:num_partition
        sys = filter(isFloat, Network[:,k])
        sys = [Int64(item)+1 for item in sys]
        partition[k]["SelfBus"] = sys
        for i in sys
            Dict_BusToZone[i] = k
        end
    end

    ## List of tie-lines
    Dict_Location_Line = Dict()
    List_BdNodes = []
    ## Generate partiton detailed information
    for k in 1:num_partition
        partition[k]["SelfLine"] = []
        partition[k]["OtherLine"] = []
        partition[k]["OtherBus"] = []
        partition[k]["BdBus"] = []
        partition[k]["SelfGen"] = []
    end
    for (i,j) in keys(Dict_Line)
        fr_zone = Dict_BusToZone[i]
        to_zone = Dict_BusToZone[j]
        if fr_zone == to_zone
            push!(partition[fr_zone]["SelfLine"], (i,j))
        else
            Dict_Location_Line[(i,j)] = [fr_zone, to_zone]
            push!(List_BdNodes, i)
            push!(List_BdNodes, j)
            push!(partition[fr_zone]["OtherLine"], (i,j))
            push!(partition[to_zone]["OtherLine"], (i,j))
            push!(partition[fr_zone]["BdBus"], i)
            push!(partition[to_zone]["BdBus"], j)
            push!(partition[fr_zone]["OtherBus"], j)
            push!(partition[to_zone]["OtherBus"], i)
        end
    end
    ## Remove duplicate
    unique!(List_BdNodes)
    for k in 1:num_partition
        unique!(partition[k]["BdBus"])
        unique!(partition[k]["OtherLine"])
        unique!(partition[k]["OtherBus"])
    end
    ## Generator
    for k in 1:num_partition
        for i in keys(Dict_Gen["GenIdx"])
            if i in partition[k]["SelfBus"]
                push!(partition[k]["SelfGen"], i)
            end
        end
    end
    ## Locations of Bdnodes
    Dict_Location = Dict(i=>[] for i in List_BdNodes)
    for i in List_BdNodes
        for k in 1:num_partition
            if (i in partition[k]["BdBus"]) || (i in partition[k]["OtherBus"])
                push!(Dict_Location[i], k)
            end
        end
    end
    return partition, Dict_Location_Line, Dict_Location
end

## Centralized Solver
function CentralizedNetworkSolver(case)
    data_network, Network = parser(case, 2)
    Dict_Nbr, List_DiLine = GetNetworkInfo(data_network)
    Dict_Bus  = data_network["Bus"]
    Dict_Line = data_network["Line"]
    Dict_Gen  = data_network["Gen"]
    B = data_network["B"]
    G = data_network["G"]
    nBus = length(Dict_Bus["BusIdx"])
    nGen = length(Dict_Gen["p_max"])

    m = Model(solver = IpoptSolver(linear_solver = "MA27"))
    @variable(m, c[i in 1:nBus], start = 1.0,
                 lowerbound = Dict_Bus["v_min"][i]^2,
                 upperbound = Dict_Bus["v_max"][i]^2)
    @variable(m, ce[List_DiLine], start = 1.0)
    @variable(m, se[List_DiLine], start = 0.0)
    @variable(m, p_g[i in 1:nGen], lowerbound = Dict_Gen["p_min"][i],
                                   upperbound = Dict_Gen["p_max"][i])
    @variable(m, q_g[i in 1:nGen], lowerbound = Dict_Gen["q_min"][i],
                                   upperbound = Dict_Gen["q_max"][i])

    for i in 1:nBus
        if i in keys(Dict_Gen["GenIdx"])
            @constraint(m, sum(p_g[k] for k in Dict_Gen["GenIdx"][i]) -
                Dict_Bus["Load_P"][i] == G[i,i]*c[i] +
                sum(G[i,j]*ce[(i,j)]-B[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
            @constraint(m, sum(q_g[k] for k in Dict_Gen["GenIdx"][i]) -
                Dict_Bus["Load_Q"][i] == -B[i,i]*c[i] +
                sum(-B[i,j]*ce[(i,j)]- G[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
        else
            @constraint(m, 0.0 - Dict_Bus["Load_P"][i] == G[i,i]*c[i] +
                sum(G[i,j]*ce[(i,j)]-B[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
            @constraint(m, 0.0 - Dict_Bus["Load_Q"][i] == -B[i,i]*c[i] +
                sum(-B[i,j]*ce[(i,j)]- G[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
        end
    end
    for (i,j) in keys(Dict_Line)
        @constraint(m, ce[(i,j)] - ce[(j,i)] == 0)
        @constraint(m, se[(i,j)] + se[(j,i)] == 0)
        @NLconstraint(m, ce[(i,j)]^2 + se[(i,j)]^2 == c[i]*c[j])
    end
    c2 = Dict_Gen["c2"]
    c1 = Dict_Gen["c1"]
    c0 = Dict_Gen["c0"]
    # obj = sum(c2[i]*p_g[i]^2 + c1[i]*p_g[i] + c0[i] for i in 1:nGen)
    obj = sum(c1[i]*p_g[i] for i in 1:nGen)
    @objective(m, Min, obj)
    return m
end

function CentralizedNetworkSolver_SOCP(case)
    data_network, Network = parser(case, 2)
    Dict_Nbr, List_DiLine = GetNetworkInfo(data_network)
    Dict_Bus  = data_network["Bus"]
    Dict_Line = data_network["Line"]
    Dict_Gen  = data_network["Gen"]
    B = data_network["B"]
    G = data_network["G"]
    nBus = length(Dict_Bus["BusIdx"])
    nGen = length(Dict_Gen["p_max"])

    m = Model(solver = IpoptSolver(linear_solver = "MA27"))
    @variable(m, c[i in 1:nBus], start = 1.0,
                 lowerbound = Dict_Bus["v_min"][i]^2,
                 upperbound = Dict_Bus["v_max"][i]^2)
    @variable(m, ce[List_DiLine], start = 1.0)
    @variable(m, se[List_DiLine], start = 0.0)
    @variable(m, p_g[i in 1:nGen], lowerbound = Dict_Gen["p_min"][i],
                                   upperbound = Dict_Gen["p_max"][i])
    @variable(m, q_g[i in 1:nGen], lowerbound = Dict_Gen["q_min"][i],
                                   upperbound = Dict_Gen["q_max"][i])

    for i in 1:nBus
        if i in keys(Dict_Gen["GenIdx"])
            @constraint(m, sum(p_g[k] for k in Dict_Gen["GenIdx"][i]) -
                Dict_Bus["Load_P"][i] == G[i,i]*c[i] +
                sum(G[i,j]*ce[(i,j)]-B[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
            @constraint(m, sum(q_g[k] for k in Dict_Gen["GenIdx"][i]) -
                Dict_Bus["Load_Q"][i] == -B[i,i]*c[i] +
                sum(-B[i,j]*ce[(i,j)]- G[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
        else
            @constraint(m, 0.0 - Dict_Bus["Load_P"][i] == G[i,i]*c[i] +
                sum(G[i,j]*ce[(i,j)]-B[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
            @constraint(m, 0.0 - Dict_Bus["Load_Q"][i] == -B[i,i]*c[i] +
                sum(-B[i,j]*ce[(i,j)]- G[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
        end
    end
    for (i,j) in keys(Dict_Line)
        @constraint(m, ce[(i,j)] - ce[(j,i)] == 0)
        @constraint(m, se[(i,j)] + se[(j,i)] == 0)
        @NLconstraint(m, ce[(i,j)]^2 + se[(i,j)]^2 <= c[i]*c[j])
    end
    c2 = Dict_Gen["c2"]
    c1 = Dict_Gen["c1"]
    c0 = Dict_Gen["c0"]
    obj = sum(c2[i]*p_g[i]^2 + c1[i]*p_g[i] + c0[i] for i in 1:nGen)
    # obj = sum(c1[i]*p_g[i] for i in 1:nGen)
    @objective(m, Min, obj)
    return m
end

## Initialization
function Initialization(partition, Dict_Location_Line, Dict_Location)
    num_partition = length(keys(partition))
    Dual_lmd = Dict(k=>Dict() for k in 1:num_partition)
    Dual_y   = Dict(k=>Dict() for k in 1:num_partition)
    Slack_z  = Dict(k=>Dict() for k in 1:num_partition)
    BlockSol = Dict(k=>Dict() for k in 1:num_partition)
    for k in 1:num_partition
        CopiedBus = union(partition[k]["OtherBus"], partition[k]["BdBus"])
        OtherLine = partition[k]["OtherLine"]
        Dual_lmd[k]["Line"] = Dict(l=>[0.0, 0.0] for l in OtherLine)
        Dual_y[k]["Line"]   = Dict(l=>[0.0, 0.0] for l in OtherLine)
        Slack_z[k]["Line"]  = Dict(l=>[0.0, 0.0] for l in OtherLine)
        BlockSol[k]["Line"] = Dict(l=>[1.0, 0.0] for l in OtherLine)

        Dual_lmd[k]["Bus"] = Dict(i=>0.0 for i in CopiedBus)
        Dual_y[k]["Bus"]   = Dict(i=>0.0 for i in CopiedBus)
        Slack_z[k]["Bus"]  = Dict(i=>0.0 for i in CopiedBus)
        BlockSol[k]["Bus"] = Dict(i=>1.0 for i in CopiedBus)
    end
    Global_copy = Dict("Line"=>Dict(l=>[0.0, 0.0] for l in keys(Dict_Location_Line)),
                       "Bus"=>Dict(i=>1.0 for i in keys(Dict_Location)))

    return Dual_lmd, Dual_y, Slack_z, BlockSol, Global_copy
end


## Construct Subproblem
@everywhere function NetworkSubproblemConstructor(partition_k::Dict, data_network::Dict, Dict_Nbr::Dict)
    Dict_Bus  = data_network["Bus"]
    Dict_Line = data_network["Line"]
    Dict_Gen  = data_network["Gen"]
    B = data_network["B"]
    G = data_network["G"]
    m = Model(solver = IpoptSolver(linear_solver="MA27", print_level=0))
    total_bus = union(partition_k["SelfBus"],  partition_k["OtherBus"])
    total_line = []
    total_gen = []
    for (i,j) in union(partition_k["SelfLine"], partition_k["OtherLine"])
        push!(total_line, (i,j))
        push!(total_line, (j,i))
    end
    for i in partition_k["SelfGen"]
        for g in data_network["Gen"]["GenIdx"][i]
            push!(total_gen, g)
        end
    end
    @variable(m, c[i in total_bus], start = 1.0,
            lowerbound = Dict_Bus["v_min"][i]^2,
            upperbound = Dict_Bus["v_max"][i]^2)
    @variable(m, ce[total_line], start = 1.0)  #max(fStart[i]-0.3,0))
    @variable(m, se[total_line], start = 0.0)
    @variable(m, cost)

    @variable(m, p_g[g in total_gen], lowerbound = data_network["Gen"]["p_min"][g],
                                      upperbound = data_network["Gen"]["p_max"][g])
    @variable(m, q_g[g in total_gen], lowerbound = data_network["Gen"]["q_min"][g],
                                      upperbound = data_network["Gen"]["q_max"][g])

    for i in partition_k["SelfBus"]
        if i in partition_k["SelfGen"]
            @constraint(m, sum(p_g[k] for k in data_network["Gen"]["GenIdx"][i]) -
                data_network["Bus"]["Load_P"][i] == G[i,i]*c[i] +
                sum(G[i,j]*ce[(i,j)]-B[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
            @constraint(m, sum(q_g[k] for k in data_network["Gen"]["GenIdx"][i]) -
                data_network["Bus"]["Load_Q"][i] == -B[i,i]*c[i] +
                sum(-B[i,j]*ce[(i,j)]- G[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
        else
            @constraint(m, 0.0 - data_network["Bus"]["Load_P"][i] == G[i,i]*c[i] +
                sum(G[i,j]*ce[(i,j)]-B[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
            @constraint(m, 0.0 - data_network["Bus"]["Load_Q"][i] == -B[i,i]*c[i] +
                sum(-B[i,j]*ce[(i,j)]- G[i,j]*se[(i,j)] for j in Dict_Nbr[i]))
        end
    end
    for (i,j) in union(partition_k["SelfLine"], partition_k["OtherLine"])
        @constraint(m, ce[(i,j)] - ce[(j,i)] == 0)
        @constraint(m, se[(i,j)] + se[(j,i)] == 0)
        @NLconstraint(m, ce[(i,j)]^2 + se[(i,j)]^2 == c[i]*c[j])
    end
    c2 = data_network["Gen"]["c2"]
    c1 = data_network["Gen"]["c1"]
    c0 = data_network["Gen"]["c0"]
    # obj_aux = sum(c2[i]*p_g[i]^2 + c1[i]*p_g[i] + c0[i] for i in total_gen)
    if length(total_gen) > 0
        obj_aux = sum(c1[i]*p_g[i] for i in total_gen)
    else
        obj_aux = 0.0
    end
    @constraint(m, obj_aux == cost)
    @objective(m, Min, cost)
    return m
end



## subproblem solvers
@everywhere function PreSolve_Network()
    global m
    solve(m)
end

@everywhere function SolveSubproblem_TL(partition_k, Global_copy, y_k, z_k, rho)
    global m
    aug_obj = m[:cost]
    for l in partition_k["OtherLine"]
        temp = temp =[m[:ce][l], m[:se][l]]-Global_copy["Line"][l]+ z_k["Line"][l]
        aug_obj =aug_obj+y_k["Line"][l]'*temp + 0.5 * rho * temp'*temp
    end
    for i in union(partition_k["BdBus"], partition_k["OtherBus"])
        temp = m[:c][i] - Global_copy["Bus"][i] + z_k["Bus"][i]
        aug_obj =aug_obj+ y_k["Bus"][i] * temp + 0.5 * rho * temp'*temp
    end
    @objective(m, Min, aug_obj)
    status = solve(m)
    gen_cost = getvalue(m[:cost])
    BlockSol_k = Dict(
        "Line"=>Dict(l=>[getvalue(m[:ce][l]), getvalue(m[:se][l])] for l in partition_k["OtherLine"]),
        "Bus"=>Dict(i=>getvalue(m[:c][i]) for i in union(partition_k["BdBus"], partition_k["OtherBus"])))
    return BlockSol_k, gen_cost
end


## Second block
function UpdateSecondBlock!(Global_copy, BlockSol, Slack_z, Dual_y, Dict_Location_Line, Dict_Location, rho)
    for l in keys(Dict_Location_Line)
        Global_copy["Line"][l] = sum(BlockSol[k]["Line"][l] + Slack_z[k]["Line"][l]
            for k in Dict_Location_Line[l])/2 +
            sum(Dual_y[k]["Line"][l] for k in Dict_Location_Line[l])/(2*rho)
    end
    for i in keys(Dict_Location)
        num_k = length(Dict_Location[i])
        Global_copy["Bus"][i] = sum(BlockSol[k]["Bus"][i] + Slack_z[k]["Bus"][i]
            for k in Dict_Location[i])/num_k +
            sum(Dual_y[k]["Bus"][i] for k in Dict_Location[i])/(num_k * rho)
    end
end

## Third block
function UpdateThirdBlock!(Slack_z, BlockSol, Global_copy, Dual_y, Dual_lmd, rho, beta)
    num_partition = length(keys(BlockSol))
    for k in 1:num_partition
        for i in keys(Slack_z[k]["Bus"])
            Slack_z[k]["Bus"][i] = (rho*(Global_copy["Bus"][i] - BlockSol[k]["Bus"][i])-
                Dual_lmd[k]["Bus"][i]-Dual_y[k]["Bus"][i])/(beta+rho)
        end
        for l in keys(Slack_z[k]["Line"])
            Slack_z[k]["Line"][l] =(rho*(Global_copy["Line"][l] - BlockSol[k]["Line"][l])-
                Dual_lmd[k]["Line"][l]-Dual_y[k]["Line"][l])/(beta+rho)
        end
    end
end
## Inner dual
function UpdateDual_y!(Dual_y, BlockSol, Global_copy,Slack_z, rho)
    num_partition = length(keys(BlockSol))
    re3_list = []
    re2_list = []
    z_list = []
    for k in 1:num_partition
        for i in keys(Slack_z[k]["Bus"])
            re3_temp = BlockSol[k]["Bus"][i]-Global_copy["Bus"][i] + Slack_z[k]["Bus"][i]
            re2_temp = BlockSol[k]["Bus"][i]-Global_copy["Bus"][i]
            Dual_y[k]["Bus"][i] =Dual_y[k]["Bus"][i]+rho*re3_temp
            push!(re3_list, re3_temp'*re3_temp)
            push!(re2_list, re2_temp'*re2_temp)
            push!(z_list, Slack_z[k]["Bus"][i]'*Slack_z[k]["Bus"][i])
        end
        for l in keys(Slack_z[k]["Line"])
            re3_temp = BlockSol[k]["Line"][l]-Global_copy["Line"][l] + Slack_z[k]["Line"][l]
            re2_temp = BlockSol[k]["Line"][l]-Global_copy["Line"][l]
            Dual_y[k]["Line"][l] =Dual_y[k]["Line"][l]+rho*re3_temp
            push!(re3_list, re3_temp'*re3_temp)
            push!(re2_list, re2_temp'*re2_temp)
            push!(z_list, Slack_z[k]["Line"][l]'*Slack_z[k]["Line"][l])
        end
    end
    re3 = sqrt(sum(re3_list))
    re2 = sqrt(sum(re2_list))
    rez = sqrt(sum(z_list))
    return re2, re3, rez
end
## Outer dual
function UpdateOuterDual!(Dual_lmd, Slack_z, beta, lmd_bd)
    num_partition = length(keys(Dual_lmd))
    for k in 1:num_partition
        for l in keys(Dual_lmd[k]["Line"])
            temp = Dual_lmd[k]["Line"][l]+beta*Slack_z[k]["Line"][l]
            Dual_lmd[k]["Line"][l] = min.(max.(temp, -lmd_bd), lmd_bd)
        end
        for i in keys(Dual_lmd[k]["Bus"])
            temp = Dual_lmd[k]["Bus"][i]+beta*Slack_z[k]["Bus"][i]
            Dual_lmd[k]["Bus"][i] = min.(max.(temp, -lmd_bd), lmd_bd)
        end
    end
end

function InitializeInner!(Dual_y, Dual_lmd, Slack_z)
    num_partition = length(keys(Dual_y))
    for k in 1:num_partition
        for l in keys(Slack_z[k]["Line"])
            Slack_z[k]["Line"][l] = [0.0, 0.0]
            Dual_y[k]["Line"][l] = -Dual_lmd[k]["Line"][l]
        end
        for i in keys(Slack_z[k]["Bus"])
            Slack_z[k]["Bus"][i] = 0.0
            Dual_y[k]["Bus"][i] = -Dual_y[k]["Bus"][i]
        end
    end
end

## ADMM_g
@everywhere function SolveSubproblem_ADMMg(partition_k, Global_copy, total_gen, y_k, z_k, rho, h)
    global m
    aug_obj = m[:cost]
    for l in partition_k["OtherLine"]
        temp = temp =[m[:ce][l], m[:se][l]]-Global_copy["Line"][l]+ z_k["Line"][l]
        aug_obj =aug_obj+y_k["Line"][l]'*temp + 0.5 * rho * temp'*temp
    end
    for i in union(partition_k["BdBus"], partition_k["OtherBus"])
        temp = m[:c][i] - Global_copy["Bus"][i] + z_k["Bus"][i]
        aug_obj =aug_obj+ y_k["Bus"][i] * temp + 0.5 * rho * temp'*temp
    end
    ## Additional prox terms
    aug_obj += 0.5 * h * sum((m[:c][i] - getvalue(m[:c][i]))^2 for
        i in union(partition_k["SelfBus"], partition_k["OtherBus"]))
    aug_obj += 0.5 * h * sum((m[:ce][l] - getvalue(m[:ce][l]))^2 for
        l in union(partition_k["SelfLine"], partition_k["OtherLine"]))
    aug_obj += 0.5 * h * sum((m[:se][l] - getvalue(m[:se][l]))^2 for
        l in union(partition_k["SelfLine"], partition_k["OtherLine"]))

    if length(total_gen) > 0
        aug_obj += 0.5 * h * sum((m[:p_g][i] - getvalue(m[:p_g][i]))^2 for i in total_gen)
        aug_obj += 0.5 * h * sum((m[:q_g][i] - getvalue(m[:q_g][i]))^2 for i in total_gen)
    end
    @objective(m, Min, aug_obj)
    status = solve(m)
    gen_cost = getvalue(m[:cost])
    BlockSol_k = Dict(
        "Line"=>Dict(l=>[getvalue(m[:ce][l]), getvalue(m[:se][l])] for l in partition_k["OtherLine"]),
        "Bus"=>Dict(i=>getvalue(m[:c][i]) for i in union(partition_k["BdBus"], partition_k["OtherBus"])))
    return BlockSol_k, gen_cost
end

function UpdateSecondBlock_ADMMg!(Global_copy, BlockSol, Slack_z, Dual_y,
    Dict_Location_Line, Dict_Location, rho, h)
    Global_copy_prev = copy(Global_copy)
    for l in keys(Dict_Location_Line)
        Global_copy["Line"][l] = (rho *sum(BlockSol[k]["Line"][l] + Slack_z[k]["Line"][l]
            for k in Dict_Location_Line[l]) +
            sum(Dual_y[k]["Line"][l] for k in Dict_Location_Line[l])+
            h*Global_copy_prev["Line"][l])/(2*rho+h)
    end
    for i in keys(Dict_Location)
        num_k = length(Dict_Location[i])
        Global_copy["Bus"][i] = (rho*sum(BlockSol[k]["Bus"][i] + Slack_z[k]["Bus"][i]
            for k in Dict_Location[i]) +
            sum(Dual_y[k]["Bus"][i] for k in Dict_Location[i]) +
            h*Global_copy_prev["Bus"][i])/(num_k * rho+h)
    end
end

function UpdateThirdBlock_ADMMg!(Slack_z, BlockSol, Global_copy, Dual_y, rho, beta, step_size)
    num_partition = length(keys(BlockSol))
    Slack_z_prev = copy(Slack_z)
    for k in 1:num_partition
        for i in keys(Slack_z[k]["Bus"])
            Slack_z[k]["Bus"][i] = Slack_z_prev[k]["Bus"][i] - step_size * (
                beta*Slack_z_prev[k]["Bus"][i] + Dual_y[k]["Bus"][i] +
                rho * (BlockSol[k]["Bus"][i]-Global_copy["Bus"][i]+Slack_z_prev[k]["Bus"][i]))
        end
        for l in keys(Slack_z[k]["Line"])
            Slack_z[k]["Line"][l] = Slack_z_prev[k]["Line"][l] - step_size * (
                beta*Slack_z_prev[k]["Line"][l] + Dual_y[k]["Line"][l] +
                rho * (BlockSol[k]["Line"][l]-Global_copy["Line"][l]+Slack_z_prev[k]["Line"][l]))
        end
    end
end

## PDD
@everywhere function SolveSubproblem_PDD(partition_k, Global_copy, y_k, rho)
    global m
    aug_obj = m[:cost]
    for l in partition_k["OtherLine"]
        temp = temp =[m[:ce][l], m[:se][l]]-Global_copy["Line"][l]
        aug_obj =aug_obj+y_k["Line"][l]'*temp + 0.5 * rho * temp'*temp
    end
    for i in union(partition_k["BdBus"], partition_k["OtherBus"])
        temp = m[:c][i] - Global_copy["Bus"][i]
        aug_obj =aug_obj+ y_k["Bus"][i] * temp + 0.5 * rho * temp'*temp
    end


    @objective(m, Min, aug_obj)
    status = solve(m)
    gen_cost = getvalue(m[:cost])
    BlockSol_k = Dict(
        "Line"=>Dict(l=>[getvalue(m[:ce][l]), getvalue(m[:se][l])] for l in partition_k["OtherLine"]),
        "Bus"=>Dict(i=>getvalue(m[:c][i]) for i in union(partition_k["BdBus"], partition_k["OtherBus"])))
    return BlockSol_k, gen_cost
end

## PDD
function SampleAndUpdateFirstBlock_PDD!(Dict_Location,
    Dict_Location_Line, BlockSol, Global_copy, Dual_y, rho)
    block_type = "Ipopt"
    block_idx = -1
    num_partition = length(keys(BlockSol))
    num_total_block = num_partition+length(keys(Dict_Location))+length(keys(Dict_Location_Line))
    idx_first_block = rand(1:num_total_block)

    if idx_first_block > num_partition && idx_first_block <= num_partition + length(keys(Dict_Location))
        i = sort(collect(keys(Dict_Location)))[idx_first_block-num_partition]
        block_type = "Bus"
        block_idx = i
        num_k = length(Dict_Location[i])
        Global_copy["Bus"][i] = sum(BlockSol[k]["Bus"][i]
            for k in Dict_Location[i])/num_k +
            sum(Dual_y[k]["Bus"][i] for k in Dict_Location[i])/(num_k * rho)
    elseif idx_first_block > num_partition + length(keys(Dict_Location))
        temp= idx_first_block -num_partition -length(keys(Dict_Location))
        l = collect(keys(Dict_Location_Line))[temp]
        block_type = "Line"
        block_idx = l
        Global_copy["Line"][l] = sum(BlockSol[k]["Line"][l]
            for k in Dict_Location_Line[l])/2 +
            sum(Dual_y[k]["Line"][l] for k in Dict_Location_Line[l])/(2*rho)
    end
    return block_type, block_idx
end

function UpdateSecondBlock_PDD!(Global_copy, BlockSol, Dual_y,
    Dict_Location_Line, Dict_Location, rho, block_idx)

    for l in keys(Dict_Location_Line)
        if l!= block_idx
            Global_copy["Line"][l] = sum(BlockSol[k]["Line"][l]
                for k in Dict_Location_Line[l])/2 +
                sum(Dual_y[k]["Line"][l] for k in Dict_Location_Line[l])/(2*rho)
        end
    end
    for i in keys(Dict_Location)
        if i != block_idx
            num_k = length(Dict_Location[i])
            Global_copy["Bus"][i] = sum(BlockSol[k]["Bus"][i]
                for k in Dict_Location[i])/num_k +
                sum(Dual_y[k]["Bus"][i] for k in Dict_Location[i])/(num_k * rho)
        end
    end
end

function Calculate_metrics!(BlockSol, Global_copy, Dual_y, gen_cost,rho)
    num_partition = length(keys(BlockSol))
    re2_list = []
    aug_obj = gen_cost
    for k in 1:num_partition
        for i in keys(BlockSol[k]["Bus"])
            re2_temp = BlockSol[k]["Bus"][i]-Global_copy["Bus"][i]
            aug_obj += Dual_y[k]["Bus"][i]'*re2_temp + 0.5 * rho * re2_temp'*re2_temp
            push!(re2_list, re2_temp'*re2_temp)

        end
        for l in keys(BlockSol[k]["Line"])
            re2_temp = BlockSol[k]["Line"][l]-Global_copy["Line"][l]
            aug_obj += Dual_y[k]["Line"][l]'*re2_temp + 0.5 * rho * re2_temp'*re2_temp
            push!(re2_list, re2_temp'*re2_temp)
        end
    end
    re2 = sqrt(sum(re2_list))
    return re2, aug_obj
end
function UpdateDual_PDD!(BlockSol, Global_copy, Dual_y,rho)
    num_partition = length(keys(BlockSol))
    for k in 1:num_partition
        for i in keys(BlockSol[k]["Bus"])
            re2_temp = BlockSol[k]["Bus"][i]-Global_copy["Bus"][i]
            Dual_y[k]["Bus"][i] = Dual_y[k]["Bus"][i] + rho*re2_temp
        end
        for l in keys(BlockSol[k]["Line"])
            re2_temp = BlockSol[k]["Line"][l]-Global_copy["Line"][l]
            Dual_y[k]["Line"][l] = Dual_y[k]["Line"][l] + rho*re2_temp
        end
    end
end


# case="case14"
# num_partition = 2
# data_network, Network = parser(case, num_partition)
# Dict_Nbr, List_DiLine = GetNetworkInfo(data_network)
# partition, Dict_Location_Line, Dict_Location = GeneratePartitionInfo(Network, num_partition, data_network)
# m = CentralizedNetworkSolver(case)
# m_SOCP = CentralizedNetworkSolver(case)
#
# solve(m)
# solve(m_SOCP)
