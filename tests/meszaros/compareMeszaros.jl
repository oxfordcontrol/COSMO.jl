# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
workspace()
include("../../src/Solver.jl")
include("ConvertProblem.jl")

using OSQP, OSSDP, Base.Test, JLD,MAT

dirPath = "./MAT_Files/"

fileNames = []
for f in filter(x -> endswith(x, ".mat"), readdir(dirPath))
    f = split(f,".")[1]
    push!(fileNames,f)
end
excludeProbs = ["BOYD1";"BOYD2";"CONT-200";"CONT-201";"CONT-300";"UBH1"]

filter!(x->!in(x,excludeProbs),fileNames)
# fileNames = fileNames[[5,6,7,8,11,12,13,21,22,23,24,25,26,28,29,30,31]]
resIter = zeros(length(fileNames),4)
resCost = zeros(length(fileNames),5)
resStatus = Array{Symbol}(length(fileNames),2)
resX = Array{Array}(length(fileNames),4)
iii = 1


fn =  "meszarosComparison.jld"



for file in fileNames

  # jump to next file if error happens
  try
    gc()
    data = matread("./MAT_FILES/$(file).mat")
    Q = data["Q"]
    A = data["A"]
    c = data["c"]
    ru = data["ru"]
    rl = data["rl"]
    lb = data["lb"]
    ub = data["ub"]


    costTrue = 0.
    try
      costTrue = JLD.load("./MAT_FILES/objVals.jld")["nameValDict"][lowercase(file)]
    catch
      costTrue = NaN
    end
    m = size(A,1)
    n = size(Q,1)

    Pa, qa, Aa, ba, K, typ = Converter.convertProblem(data)
    #println("Problem type:$(typ == 1 ? "Ax=b" : "l<=Ax<=u")")
    settings = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=1500,verbose=false,scaling = 0,checkTermination = 1,eps_abs = 1e-3,eps_rel=1e-4,timelimit=60)
    resOSSDP_unscaled,nothing = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    print(".")
    settings = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=1500,verbose=false,scaling = 10,checkTermination = 1,eps_abs = 1e-3,eps_rel=1e-4,timelimit=60)
    resOSSDP_scaled,nothing = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    print(".")


    # solve with OSQP
    Aa2 = [A;eye(n)]
    P2 = Q
    q2 = c
    if norm(rl-ru,Inf) < 1e-4
      u = [rl;ub]
      l = [rl;lb]
    else
      l = [rl;lb]
      u = [ru;ub]
    end

    m2u = OSQP.Model()
    OSQP.setup!(m2u; P=P2, q=q2[:], A=Aa2, l=l[:], u=u[:],scaling=0,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
    resOSQP_unscaled = OSQP.solve!(m2u)
    print(".")

    m2s = OSQP.Model()
    OSQP.setup!(m2s; P=P2, q=q2[:], A=Aa2, l=l[:], u=u[:],scaling=10,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
    resOSQP_scaled = OSQP.solve!(m2s)
    print(".")

    resCost[iii,:] = [resOSSDP_unscaled.cost resOSSDP_scaled.cost resOSQP_unscaled.info.obj_val resOSQP_scaled.info.obj_val costTrue]
    resIter[iii,:] = [resOSSDP_unscaled.iter resOSSDP_scaled.iter resOSQP_unscaled.info.iter resOSQP_scaled.info.iter]
    resStatus[iii,:] = [resOSSDP_unscaled.status resOSSDP_scaled.status]
    resX[iii,1] = resOSSDP_unscaled.x[1:n]
    resX[iii,2] = resOSSDP_scaled.x[1:n]
    resX[iii,3] = resOSQP_unscaled.x
    resX[iii,4] = resOSQP_scaled.x

    println("$(iii)/$(length(fileNames)) $(file) completed! (unscaled status: $(resOSSDP_unscaled.status), scaled status: $(resOSSDP_scaled.status))")
    iii +=1
    JLD.save(fn, "resCost", resCost, "resIter",resIter,"resX",resX,"fileNames",fileNames)

  catch
    println("An error happened with file $(file).")
    continue
  end
end


println(">>> Test Data successfully saved in $(fn).\n")
