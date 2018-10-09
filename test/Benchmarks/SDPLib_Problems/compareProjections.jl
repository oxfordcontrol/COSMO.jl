# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
include("../../../src/QOCS.jl")

using Main.QOCS, Test, JLD2
using SparseArrays
using JuMP, SCS

dirPath = "../sdplib/"

fileNames = String[]
matrixDims = Int[]
for f in filter(x -> endswith(x, ".jld"), readdir(dirPath))
    f = split(f,".")[1]
    @load dirPath*f*".jld" n
    push!(fileNames,f)
    push!(matrixDims, n)
end


# sort filenames by number of nnz (stored in problemData[:,4])
sortedInd = sort!(collect(1:1:length(fileNames)), by=i->matrixDims[i])
fileNames = fileNames[sortedInd]

# filter some problems by name
# excludeProbs = ["BOYD1";"BOYD2";"CONT-200";"CONT-201";"CONT-300";"UBH1"]
# filter!(x->!in(x,excludeProbs),fileNames)

# to begin with only look at first nn problems
pushfirst!(fileNames, fileNames[1])
fileNames = fileNames[1:30]

outputFilename = "sdplibProjectionComparison.jld"
count = 1;
for file in fileNames
  # jump to next file if error happens
  @show file
  @load  dirPath*file*".jld" m n nblocks blockvec c F optVal

  P = spzeros(m,m)
  q = Array(c)  # Make dense
  A = spzeros(n^2,m)
  for iii = 1:m
    A[:,iii] = -vec(F[iii+1])
  end
  b = Array(-vec(F[1]))  # Make dense

  constraint1 = QOCS.Constraint(-A,b,QOCS.PositiveSemidefiniteCone())
  # settings = QOCS.Settings(sigma=1e-6,alpha=1.6,max_iter=5000,verbose=true,check_termination=1,eps_abs = 1e-6, eps_rel = 1e-6)
  settings = QOCS.Settings(adaptive_rho=false,verbose=true, use_lanczos=false, max_iter=5000)
  model = QOCS.Model()
  assemble!(model,P,q,(constraint1))
  @time res = QOCS.optimize!(model,settings);
  @show res
  
  mosek_model = JuMP.Model()
  setsolver(mosek_model, SCSSolver(verbose=false))
  @variable(mosek_model, x[1:m])
  @variable(mosek_model, S[1:n,1:n],SDP)
  @objective(mosek_model, Min, q'*x)
  s = vec(S)
  @constraint(mosek_model, A*x+s .== b)
  status = JuMP.solve(mosek_model)

  # println("$(count)/$(length(fileNames)) $(file) completed! (, scaled status: $(resOSSDP_scaled.status))")
  global count = count + 1
  #JLD.save(outputFilename, "resCost", resCost, "resIter",resIter,"resX",resX,"fileNames",fileNames)
end
