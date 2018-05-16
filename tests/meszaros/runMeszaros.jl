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
readmeInfo = JLD.load("./bart_meszaros_data/objVals.jld")
problemData = readmeInfo["problemData"]
sortedInd = sort!(collect(1:1:length(fileNames)), by=i->problemData[i,4])
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

for file in fileNames
  # jump to next file if error happens
  # try
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

    Pa, qa, r, Aa, ba, K = Converter.convertProblem(data)
    settings = OSSDPSettings(adaptive_rho=true, max_iter=4000)
    print("Running QOCS:")
    @time res, nothing = OSSDP.solve(Pa,qa,Aa,ba,K,settings)

    m_ = OSQP.Model()
    OSQP.setup!(m_; P=P, q=q[:], A=A, l=l[:], u=u[:], verbose=false)
    print("Running OSQP:")
    @time resOSQP = OSQP.solve!(m_)
    print(".")

    # add the constant term to the solution
    resCost[iii,:] = [res.cost costTrue]
    resIter[iii,:] = [res.iter]
    resStatus = Array{Symbol}(length(fileNames),1)

    println("Diff OSQP-QOCS (scaled): $(100*(res.cost - resOSQP.info.obj_val)/resOSQP.info.obj_val)%")
    println("Iter OSQP-QOCS (scaled): $(100*(res.iter - resOSQP.info.iter)/resOSQP.info.iter)%")

    println("$(iii)/$(length(fileNames)) $(file) completed! (status: $(res.status))")
    iii +=1
    JLD.save(fn, "resCost", resCost, "resIter",resIter, "fileNames", fileNames)

  #=
  catch
    println("An error happened with file $(file).")
    continue
  end
  =#
end


println(">>> Test Data successfully saved in $(fn).\n")
