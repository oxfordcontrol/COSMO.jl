# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
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
readmeInfo = JLD.load("./bart_meszaros_data/objVals.jld")
problemData = readmeInfo["problemData"]
nnzs = problemData[:,4] + problemData[:,6]
sortedInd = sort!(collect(1:1:length(fileNames)), by=i->nnzs[i])
fileNames = fileNames[sortedInd]

# filter some problems by name
#=
excludeProbs = ["BOYD1";"BOYD2";"CONT-200";"CONT-201";"CONT-300";"UBH1"]
filter!(x->!in(x,excludeProbs),fileNames)
=#

# to begin with only look at first nn problems
fileNames = fileNames[120:end]

resCost = zeros(length(fileNames),2)
resIter = zeros(length(fileNames),1)
resStatus = zeros(length(fileNames),1)
# resTime = zeros(length(fileNames),1)
iii = 1



timestamp = Dates.format(now(), "yyddmm_HH-MM")
fn = timestamp * "meszarosComparison.jld"

for file in ["DTOC3", "DTOC3"] # ["QSCTAP2", "QSCTAP2"] # fileNames
  # jump to next file if error happens
  println("----------------------------------------")
  print(file)
  flush(STDOUT)
  # local data, Pa, Aa, res, tt
  # let data, Pa, Aa, res, tt
    data = matread("$(dirPath)"*"$(file).mat")
    r = data["r"]
    costTrue = 0.
    try
      costTrue = readmeInfo["nameValDict"][lowercase(file)]
    catch
      costTrue = NaN
    end

    Pa, qa, r, Aa, ba, K = Converter.convertProblem(data)
    qa = full(qa[:, 1])
    ba = full(ba)
    println("  |  nnz: $(nnz(Pa) + nnz(Aa))")
    println("----------------------------------------")

    settings = OSSDPSettings(adaptive_rho=true, max_iter=5000, verbose=false, checkTermination=100, scaling=10)
    # print("Running QOCS:")
    res, tt = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    # Profile.clear()
    # @profile res, tt = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    # open(Profile.print, "profile.txt", "w")

    m_ = OSQP.Model()
    OSQP.setup!(m_; P=data["P"], q=data["q"][:], A=data["A"], l=data["l"][:], u=data["u"][:], verbose=false, max_iter=5000, check_termination=100, adaptive_rho_interval=1)
    # print("Running OSQP:")
    resOSQP = OSQP.solve!(m_)

    println("OSQP avg iter time: $(resOSQP.info.run_time/resOSQP.info.iter)")
    println("Diff Obj (QOCS-OSQP)/OSQP: $(100*(res.cost - resOSQP.info.obj_val)/resOSQP.info.obj_val)%")
    println("QOCS iter: $(res.iter) [ratio QOCS/OSQP $(res.iter/resOSQP.info.iter)]")

    println("$(iii)/$(length(fileNames)) $(file) completed! (status QOCS: $(res.status), status OSQP: $(resOSQP.info.status).)")
    iii +=1
  # end
end
