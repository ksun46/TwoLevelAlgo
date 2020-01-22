## Robust tensor PCA problem with l_1 norm regularization
## ADMM-g and Two-level Algorithm

include("../warmup.jl")
import LinearAlgebra:norm

## Generate data with R = Rcp + ceil(0.2*Rcp)
function GenerateData(Rcp)
    ## Construct Problem Data
    ## Construct Z as outer product of Gaussian vectors
    Z = sum(ttt(ttt(rand(I1), rand(I2)), rand(I3)) for i in 1:Rcp)
    ## Construct E as a sparse tensor with sparsity 0.001*I1*I2*I3
    E = zeros(I1,I2,I3)
    for t in 1: trunc(Int64,0.001*I1*I2*I3)
        i = rand(1:I1)
        j = rand(1:I2)
        k = rand(1:I3)
        E[i,j,k] = randn()
    end
    ## Construct BB as 0.001 * Gaussian entries
    BB = randn(I1,I2,I3) * 0.001

    ## Mode-1 unfolding of T
    T = Z + E + BB
    T1 = tenmat(T,1)

    ## Initial guess
    R = trunc(Int64,Rcp + ceil(0.2 * Rcp)) ## Initial guess of CP rank
    A = randn(I1, R)
    B = randn(I2, R)
    C = randn(I3, R)
    Z1 = tenmat(Z,1)
    E1 = tenmat(E,1)
    BB1 =tenmat(BB,1)
    Y1 = tenmat(zeros(I1,I2,I3),1) ## mode-1 unfolfing of ADMM dual variable Y
    return Z, E, BB, T, T1, R, A, B, C, Z1,E1,BB1, Y1
end

## Khatri-Rao product: column-wise Kronecker Product
function KR(A::Array{Float64, 2}, B::Array{Float64,2})
    N_1 = size(A)[1]
    N_2 = size(B)[1]
    N = size(A)[2]
    Z = zeros(N_1*N_2, N)
    for i in 1:N
        Z[:, i] = kron(A[:,i],B[:,i])
    end
    return Z
end

## Soft Shrinkage Operator
function SoftShrinkage(Input::Array{Float64, 2}, a::Float64)
    n = size(Input)[1]
    m = size(Input)[2]
    T = zeros(n,m)
    for i in 1:n, j in 1:m
        val = Input[i,j]
        if abs(val) > a
            T[i,j] = sign(val) * (abs(val)-a)
        end
    end
    return T
end

## ADMM_g
function ADMM_g(A, B, C, E1, Z,Z1, BB1, T1, Y1, R)
    global I1, I2, I3, alpha, alphaN, res_tol,  maxIter
    A = deepcopy(A)
    B = deepcopy(B)
    C = deepcopy(C)
    E1 = deepcopy(E1)
    Z1 = deepcopy(Z1)
    Z_Start = deepcopy(Z1)
    BB1 = deepcopy(BB1)
    T1 = deepcopy(T1)
    Y1 = deepcopy(Y1)

    ## Parameters for ADMM-g according to Jiang et al. 2019
    rho = 4.0             ## ADMM penalty
    stepSize = 1/rho    ## Last block stpe size
    m = rho/3           ## proximal coefficient
    iterCount = 1
    reList = Float64[]
    errList = Float64[]
    objList = Float64[]
    while iterCount <= maxIter
        # Primal blocks update
        A = (tenmat(Z,1)*KR(C,B)+0.5*m*A) * inv((C'*C).*(B'*B)+0.5*m*eye(R))
        B = (tenmat(Z,2)*KR(C,A)+0.5*m*B) * inv((C'*C).*(A'*A)+0.5*m*eye(R))
        C = (tenmat(Z,3)*KR(B,A)+0.5*m*C) * inv((B'*B).*(A'*A)+0.5*m*eye(R))
        E1 = SoftShrinkage((rho/(rho+m))*(T1+(1/rho)*Y1-BB1-Z1)+ m/(rho+m)*E1, alpha/(rho+m))
        Z1 = (1/(2+2*m+rho)) * (2*A*(KR(C,B))' + 2*m*Z1 + Y1-rho*(E1+BB1-T1))
        BB1 = BB1 - stepSize*(2*alphaN*BB1-Y1+rho*(E1+Z1+BB1-T1))
        # BB1 = (Y1-rho*(Z1+E1-T1))/(rho+2*alphaN)

        # Dual update
        Y1 = Y1 - rho*(Z1+E1+BB1-T1)
        # Obj
        obj = norm(tenmat(Z,1) - A*KR(C,B)', 2)^2 + alpha*norm(E1,1) + alphaN*norm(BB1,2)^2
        push!(objList, obj)
        re = norm(E1+Z1+BB1-T1)
        push!(reList, re)
        error = norm(Z1-Z_Start)/norm(Z_Start)
        push!(errList, error)
        println("Iteration $iterCount: current residual is $re; error is $error")
        # if re <= res_tol
        #     println("ADMM_G terminate successfully!")
        #     break
        # end
        iterCount += 1
    end
    return Z1, reList, errList
end

## Two-level Algorithm with ADMM-g in the inner level
function Two_Level_g(A, B, C, E1, Z, Z1, BB1, T1, Y1, R)
    global I1, I2, I3, alpha, alphaN, res_tol, maxIter
    A = deepcopy(A)
    B = deepcopy(B)
    C = deepcopy(C)
    Z_Start = deepcopy(Z1)
    E1 = deepcopy(E1)
    Z1 = deepcopy(Z1)
    BB1 =deepcopy(BB1)
    T1 = deepcopy(T1)
    Y1 = deepcopy(Y1)
    ## Parameters for ADMM-g
    beta = 2
    rho = 4
    stepSize = 1/rho
    m = 0.5 * rho
    iterCount = 1
    almCount = 1
    ## Outer level parameters
    Lmd1 = Y1*0.0
    S1 = Y1*0.0
    reList = Float64[]
    errList = Float64[]
    flag_term = false
    while iterCount <= maxIter
        A = (tenmat(Z,1)*KR(C,B)+0.5*m*A) * inv((C'*C).*(B'*B)+0.5*m*eye(R))
        B = (tenmat(Z,2)*KR(C,A)+0.5*m*B) * inv((C'*C).*(A'*A)+0.5*m*eye(R))
        C = (tenmat(Z,3)*KR(B,A)+0.5*m*C) * inv((B'*B).*(A'*A)+0.5*m*eye(R))
        E1 = SoftShrinkage((rho/(rho+m))*(T1+(1/rho)*Y1-BB1-Z1-S1)+ m/(rho+m)*E1, alpha/(rho+m))
        Z1 = (1/(2+2*m+rho)) * (2*A*(KR(C,B))' + 2*m*Z1 + Y1-rho*(E1+BB1+S1-T1))
        # BB1 = (Y1-rho*(Z1+E1+S1-T1))/(rho+2*alphaN)
        BB1 = (Y1+m*BB1-rho*(Z1+E1+S1-T1))/(2*alphaN+rho+m)
        S1 = S1 - stepSize*(-Y1+Lmd1+rho*(Z1+E1+BB1+S1-T1)+beta*S1)
        Y1 = Y1 - rho*(Z1+E1+BB1+S1-T1)

        re3 = norm(E1+Z1+BB1+S1-T1)
        re2 = norm(E1+Z1+BB1-T1)
        error = norm(Z1-Z_Start)/norm(Z_Start)
        push!(reList, re2)
        push!(errList, error)
        println("Iteration $iterCount: current inner residual is $re3; error is $error")
        if re3 <= max(1e-3/almCount, res_tol)
            println("Inner ADMM_G terminate successfully!")
            if re2 <= res_tol
                println("Outer ADMM terminates successfully at iteration $almCount")
                flag_term = true
                # break
            end
            Lmd1 = Lmd1 + beta * S1
            if beta * 1.5 < 1e6
                beta = 1.5 * beta
                rho = 3*beta
            end
            stepSize = 1/(rho)
            almCount += 1
        end
        iterCount += 1
    end
    println("Total number of inner iterations: $iterCount")
    println("Total number of outer iterations: $almCount")
    println("      Penalty rho at termination: $rho")
    return Z1, reList, errList
end

## Dimension of the third order tensor and CP rank
I1 = 30
I2 = 50
I3 = 70
alpha = 2/max(sqrt(I1), sqrt(I2), sqrt(I3)) ## l_1 coefficient
alphaN = 1.0 ## F-norm coefficent
res_tol = 1e-5
maxIter = 2000
solDict = Dict()
RcpList = [8, 10, 20, 40]
num_test = 10
for Rcp in RcpList
    tempDict = Dict()
    err=[]
    resADMM = []
    resTL = []
    errADMM = []
    errTL = []
    for k in 1:num_test
        Z, E, BB, T, T1, R, A, B, C, Z1,E1,BB1, Y1 = GenerateData(Rcp)
        Z_one, reList1, errList1 = ADMM_g(A, B, C, E1, Z,Z1, BB1, T1, Y1, R)
        Z_two, reList2, errList2 = Two_Level_g(A, B, C, E1, Z,Z1, BB1, T1, Y1, R)
        error1 =norm(Z_one-Z1)/norm(Z1)
        error2 =norm(Z_two-Z1)/norm(Z1)
        push!(err, (error1, error2))
        push!(resADMM, reList1)
        push!(resTL, reList2)
        push!(errADMM, errList1)
        push!(errTL, errList2)
    end
    tempDict["errList"] = err
    tempDict["resADMM"] = resADMM
    tempDict["resTL"] = resTL
    tempDict["errADMM"] = errADMM
    tempDict["errTL"] = errTL
    solDict[Rcp] = tempDict
end

plotList = []
for Rcp in RcpList
    admmGeo = Float64[0.0 for i in 1:maxIter]
    twolevelGeo = Float64[0.0 for i in 1:maxIter]
    for i in 1:maxIter
        val_admm = 1.0
        val_tl = 1.0
        for k in 1:num_test
            val_admm = val_admm * solDict[Rcp]["resADMM"][k][i]
            val_tl = val_tl * solDict[Rcp]["resTL"][k][i]
        end
        admmGeo[i] = val_admm^(1/num_test)
        twolevelGeo[i] = val_tl^(1/num_test)
    end
    reList1_log = Float64[log10(admmGeo[i]) for i in 1:maxIter]
    reList2_log = Float64[log10(twolevelGeo[i]) for i in 1:maxIter]
    p = plot(reList1_log, label = "ADMM-g", linestyle = :dot, xlabel = "Iteration (CP Rank=$Rcp)", ylabel = "Infeasibility")
    plot!(reList2_log, label = "Two-level")
    push!(plotList, p)
end
p1 = plotList[1]
p2 = plotList[2]
p3 = plotList[3]
p4 = plotList[4]
p_combine = plot(p1, p2, p3, p4, layout = (2,2))
plot!(legendfont=Plots.font(5))
plot!(legend=:best)
savefig(p_combine, "combined_pca_res.pdf")

errPlotList = []
for Rcp in RcpList
    admmGeo = Float64[0.0 for i in 1:maxIter]
    twolevelGeo = Float64[0.0 for i in 1:maxIter]
    for i in 1:maxIter
        val_admm = 1.0
        val_tl = 1.0
        for k in 1:num_test
            val_admm = val_admm * solDict[Rcp]["errADMM"][k][i]
            val_tl = val_tl * solDict[Rcp]["errTL"][k][i]
        end
        admmGeo[i] = val_admm^(1/num_test)
        twolevelGeo[i] = val_tl^(1/num_test)
    end
    p = plot(admmGeo, label = "ADMM-g",linestyle = :dot, xlabel = "Iteration (CP Rank=$Rcp)", ylabel = "Rel. Error")
    plot!(twolevelGeo, label = "Two-level")
    push!(errPlotList, p)
end

e1 = errPlotList[1]
e2 = errPlotList[2]
e3 = errPlotList[3]
e4 = errPlotList[4]
e_combine = plot(e1, e2, e3, e4, layout = (2,2))
plot!(legendfont=Plots.font(5))
plot!(legend=:best)
savefig(e_combine, "combined_pca_err.pdf")
