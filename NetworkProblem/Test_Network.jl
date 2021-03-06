include("Twolevel_Network.jl")
include("PDD_Network.jl")
include("ADMMg_Network.jl")
using Dates

function WriteNetworkSol(solDict, List_network, List_partition)
    now = Dates.format(Dates.now(), "mmdd-HHMM")
    output_dir = "NetworkResult$now.csv"
    open(output_dir, "w") do f
        write(f, "case, Central_Obj, SOCP_LB, Method, Outer, Inner, Res, Obj, Gap, Time\n")

        for case in List_network
            for k in List_partition
                case_k = "$(case)-$(k)"
                central_obj = solDict[case]["central"]
                central_socp= solDict[case]["socp"]

                admm_inner = solDict[case][k]["ADMMg"]["Inner"]
                admm_res = solDict[case][k]["ADMMg"][ "Res"]
                admm_obj = solDict[case][k]["ADMMg"]["Cost"]
                admm_time = solDict[case][k]["ADMMg"]["Time"]
                admm_gap = (admm_obj-central_socp)/admm_obj
                write(f, "$case_k, $central_obj, $central_socp, ADMMg, , $admm_inner, $admm_res, $admm_obj, $admm_gap, $admm_time\n")

                pdd_inner = solDict[case][k]["PDD"]["Inner"]
                pdd_outer = solDict[case][k]["PDD"]["Outer"]
                pdd_res = solDict[case][k]["PDD"][ "Res"]
                pdd_obj = solDict[case][k]["PDD"]["Cost"]
                pdd_time = solDict[case][k]["PDD"]["Time"]
                pdd_gap = (pdd_obj - central_socp)/pdd_obj
                write(f, " , , , PDD, $pdd_outer, $pdd_inner, $pdd_res, $pdd_obj,$pdd_gap, $pdd_time\n")

                tl_inner = solDict[case][k]["Proposed"]["Inner"]
                tl_outer = solDict[case][k]["Proposed"]["Outer"]
                tl_res = solDict[case][k]["Proposed"][ "Res"]
                tl_obj = solDict[case][k]["Proposed"]["Cost"]
                tl_time = solDict[case][k]["Proposed"]["Time"]
                tl_gap = (tl_obj-central_socp)/tl_obj
                write(f, " , , , Proposed, $tl_outer, $tl_inner, $tl_res, $tl_obj,$tl_gap, $tl_time\n")
            end
        end
    end
end
# test instances
List_network = ["case14", "case118", "case300", "case1354"]
List_partition = [2, 3, 4]

solDict = Dict()
for case in List_network
    caseDict = Dict()
    ## run centralized Ipopt on nonconvex problem
    m_central = CentralizedNetworkSolver(case)
    solve(m_central)
    ## run centralized Ipopt on convex relaxation to obtain lower bound
    m_socp = CentralizedNetworkSolver_SOCP(case)
    solve(m_socp)
    ## record centralized solution
    central_obj = getobjectivevalue(m_central)
    socp_obj = getobjectivevalue(m_socp)
    caseDict["central"] = central_obj
    caseDict["socp"] = socp_obj

    for k in List_partition ## for different number of partitions of the network
        caseDict[k] = Dict()
        ## run ADMM-g
        caseDict[k]["ADMMg"] = ADMMg_Network(case, k)
        ## run PDD
        caseDict[k]["PDD"] = PDD_Network(case, k)
        ## run Two-level
        caseDict[k]["Proposed"] = Twolevel_Network(case, k)
    end
    solDict[case] = caseDict
end

WriteNetworkSol(solDict, List_network, List_partition)
