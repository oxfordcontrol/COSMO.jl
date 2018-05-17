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
fileNames = fileNames[1:end]

resCost = zeros(length(fileNames),2)
resIter = zeros(length(fileNames),1)
resStatus = zeros(length(fileNames),1)
# resTime = zeros(length(fileNames),1)
iii = 1



timestamp = Dates.format(now(), "yyddmm_HH-MM")
fn = timestamp * "meszarosComparison.jld"

for file in ["QSCTAP2"]
  # jump to next file if error happens
  println("----------------------------------------")
  print(file)
  flush(STDOUT)
  local data, Pa, Aa, res, tt
  let data, Pa, Aa, res, tt
    data = matread("$(dirPath)"*"$(file).mat")
    r = data["r"]
    costTrue = 0.
    try
      costTrue = readmeInfo["nameValDict"][lowercase(file)]
    catch
      costTrue = NaN
    end

    Pa, qa, r, Aa, ba, K = Converter.convertProblem(data)
    # qa = full(qa[:, 1])
    # ba = full(ba)
    println("  |  nnz: $(nnz(Pa) + nnz(Aa))")
    println("----------------------------------------")

    settings = OSSDPSettings(adaptive_rho=true, max_iter=4000, verbose=false, checkTermination=20)
    print("Running QOCS:")
    @time res, tt = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    # Profile.clear()
    # @profile res, tt = OSSDP.solve(Pa,qa,Aa,ba,K,settings)
    # open(Profile.print, "profile.txt", "w")

    m_ = OSQP.Model()
    OSQP.setup!(m_; P=data["P"], q=data["q"][:], A=data["A"], l=data["l"][:], u=data["u"][:], verbose=false)
    print("Running OSQP:")
    @time resOSQP = OSQP.solve!(m_)

    # add the constant term to the solution
    resCost[iii,:] = [res.cost costTrue]
    resIter[iii,:] = [res.iter]
    resStatus = Array{Symbol}(length(fileNames),1)

    println("Diff OSQP-QOCS (scaled): $(100*(res.cost - resOSQP.info.obj_val)/resOSQP.info.obj_val)%")
    println("Iter OSQP-QOCS (scaled): $(100*(res.iter - resOSQP.info.iter)/resOSQP.info.iter)%")

    println("$(iii)/$(length(fileNames)) $(file) completed! (status: $(res.status))")
    # JLD.save(fn, "resCost", resCost, "resIter",resIter, "fileNames", fileNames)
    iii +=1
  end
end


println(">>> Test Data successfully saved in $(fn).\n")
