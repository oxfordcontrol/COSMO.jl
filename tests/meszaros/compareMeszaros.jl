# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
workspace()
include("../../src/Solver.jl")
include("ConvertProblem.jl")

using OSQP, OSSDP, Base.Test, JLD,MAT

dirPath = "./bart_meszaros_data/"

fileNames = []
for f in filter(x -> endswith(x, ".mat"), readdir(dirPath))
    f = split(f,".")[1]
    push!(fileNames,f)
end


# sort filenames by number of nnz (stored in problemData[:,4])
readmeInfo = JLD.load("./MAT_FILES/objVals.jld")
problemData = readmeInfo["problemData"]
sortedInd = sort!(collect(1:1:length(fileNames)), by=i->problemData[i,4])
fileNames = fileNames[sortedInd]

# filter some problems by name
excludeProbs = ["BOYD1";"BOYD2";"CONT-200";"CONT-201";"CONT-300";"UBH1"]
filter!(x->!in(x,excludeProbs),fileNames)

# to begin with only look at first nn problems
fileNames = fileNames[1:10]

resIter = zeros(length(fileNames),4)
resTime = zeros(length(fileNames),4)
resCost = zeros(length(fileNames),5)
resStatus = Array{Symbol}(length(fileNames),2)
resX = Array{Array}(length(fileNames),4)
iii = 1



timestamp = Dates.format(now(), "yyddmm_HH-MM")
fn = timestamp * "meszarosComparison.jld"

for file in fileNames

  # jump to next file if error happens
  #try
    gc()
    data = matread("$(dirPath)"*"$(file).mat")
    P = data["P"]
    A = data["A"]
    q = data["q"]
    u = data["u"]
    l = data["l"]
    r = data["r"]

    costTrue = 0.
    try
      costTrue = readmeInfo["nameValDict"][lowercase(file)]
    catch
      costTrue = NaN
    end
    m = size(A,1)
    n = size(P,1)

    Pa, qa, Aa, ba, K = Converter.convertProblem(data)
    settings = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=1500,verbose=false,scaling = 0,checkTermination = 1,eps_abs = 1e-3,eps_rel=1e-3,timelimit=60)
    resOSSDP_unscaled,nothing = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    print(".")
    settings = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=1500,verbose=false,scaling = 10,checkTermination = 1,eps_abs = 1e-3,eps_rel=1e-3,timelimit=60)
    resOSSDP_scaled,nothing = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    print(".")


    # solve with OSQP
    m2u = OSQP.Model()
    OSQP.setup!(m2u; P=P, q=q[:], A=A, l=l[:], u=u[:],scaling=0,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
    resOSQP_unscaled = OSQP.solve!(m2u)
    print(".")

    m2s = OSQP.Model()
    OSQP.setup!(m2s; P=P, q=q[:], A=A, l=l[:], u=u[:],scaling=10,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
    resOSQP_scaled = OSQP.solve!(m2s)
    print(".")

    # add the constant term to the solution
    resCost[iii,:] = [resOSSDP_unscaled.cost+r resOSSDP_scaled.cost+r resOSQP_unscaled.info.obj_val+r resOSQP_scaled.info.obj_val+r costTrue]
    resIter[iii,:] = [resOSSDP_unscaled.iter resOSSDP_scaled.iter resOSQP_unscaled.info.iter resOSQP_scaled.info.iter]
    resStatus[iii,:] = [resOSSDP_unscaled.status resOSSDP_scaled.status]
    resX[iii,1] = resOSSDP_unscaled.x[1:n]
    resX[iii,2] = resOSSDP_scaled.x[1:n]
    resX[iii,3] = resOSQP_unscaled.x
    resX[iii,4] = resOSQP_scaled.x

    println("$(iii)/$(length(fileNames)) $(file) completed! (unscaled status: $(resOSSDP_unscaled.status), scaled status: $(resOSSDP_scaled.status))")
    iii +=1
    #JLD.save(fn, "resCost", resCost, "resIter",resIter,"resX",resX,"fileNames",fileNames)

  # catch
  #   println("An error happened with file $(file).")
  #   continue
  # end
end


println(">>> Test Data successfully saved in $(fn).\n")
