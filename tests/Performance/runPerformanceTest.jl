# Run a few test problems and determine solver time and print results to a text file (for comparison)

workspace()
include("../../src/Solver.jl")
include("../solverComparison/Compare.jl")

using OSSDP, Base.Test, Compare, JLD

BRANCH = "testMergeRequest-33"
COMMIT = ""
CONFIGURATION = "Obvious improvements only"
timestamp = Dates.format(now(), "yyddmm_HH-MM-SS")
rng = MersenneTwister(12345)

fn = "results-$(BRANCH)-$(timestamp).txt"
open(fn, "w") do f
  write(f, "BRANCH: $(BRANCH)\nCOMMIT: $(COMMIT)\nCONFIGURATION: $(CONFIGURATION)\nTIME STAMP: $(timestamp)\n\n")
end

SAVE_ALWAYS = false
maxIter = 2000
# ---------------------------------------
# RUN SOME SDP BENCHMARK PROBLEMS
# ---------------------------------------
println("> RUN SDP Benchmark Problems")
# find available problem types in SDP_Benchmark_Problems folder
probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Julia/"
existingFolders = readdir(probPath)
problemTypes = []
for f in filter(x -> !startswith(x, "."), readdir(probPath))
    f = split(f,".")[1]
    push!(problemTypes,String(f))
end
filter!(x->!in(x,["MIQO";"Lovasz";"SmallestCircle"]),problemTypes)
#filter!(x->in(x,["SDPQuad"]),problemTypes)
println(">>> $(length(problemTypes)) Problem Type(s) detected!")

kkk = 1
MeanData = zeros(1000,4)

# run tests for each problem type
for pType in problemTypes
  println(">>> ProblemType: $(pType)")
  dirPath = probPath*pType*"/"

  # find all file names in problem type folder
  problems = []
  for f in filter(x -> endswith(x, ".jld"), readdir(dirPath))
      f = split(f,".")[1]
      push!(problems,String(f))
  end
  nn = length(problems)
  # take only 5 random problems per problemtype
  rInd = randperm(rng,10)[1:3]
  problems = problems[rInd]

  problemType = JLD.load("$(dirPath)"*"$(problems[1]).jld","problemType")
  sr1 = SolverResult(nn, problemType,"QOCS - Un - noRho",timestamp,0,false,false)
  sr2 = SolverResult(nn, problemType,"QOCS - Avg - noRho",timestamp,0,true,false)
  sr3 = SolverResult(nn, problemType,"QOCS - Un",timestamp,0,false,true)
  sr4 = SolverResult(nn, problemType,"QOCS - Avg",timestamp,0,true,true)

  resData = [sr1;sr2;sr3;sr4]
  # loop over random problems in problem type folder
  for iii =1:1:length(problems)
    iii == length(problems) && (SAVE_ALWAYS = true)
    gc()
    problem = problems[iii]
    data = JLD.load("$(dirPath)"*"$(problem).jld")
    P = data["P"]
    q = vec(data["q"])
    r = data["r"]
    A = data["A"]
    b = vec(data["b"])
    m = data["m"]
    n = data["n"]
    Kf = data["Kf"]
    Kl = data["Kl"]
    Kq = data["Kq"]
    Ks = data["Ks"]

    objTrue = data["objTrue"]
    problemName = data["problemName"]


    pDims = [size(A,1);size(A,2);nnz(A)]
    # update the true value for the QOCP solver
    setUnNonAdaptive =   OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 0 ,scaleFunc=1,adaptive_rho=false)
    setMeanNonAdaptive = OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 10,scaleFunc=2,adaptive_rho=false)
    setUnAdaptive =   OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 0 ,scaleFunc=1,adaptive_rho=true)
    setMeanAdaptive = OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 10,scaleFunc=2,adaptive_rho=true)

    # define cone membership
    K = Cone(Kf,Kl,Kq,Ks)


    # Solve with OSSDP
    println(problemName)
    res1,ws1 = OSSDP.solve(P,q,A,b,K,setUnNonAdaptive);
    print("$(sr1.solverName): Cost: $(res1.cost), AvgTime/Iter: $(res1.avgIterTime) ")
    res2,ws2 = OSSDP.solve(P,q,A,b,K,setMeanNonAdaptive);
    print("$(sr2.solverName): Cost: $(res2.cost), AvgTime/Iter: $(res2.avgIterTime) ")
    res3,ws3 = OSSDP.solve(P,q,A,b,K,setUnAdaptive);
    print("$(sr3.solverName): Cost: $(res3.cost), AvgTime/Iter: $(res3.avgIterTime) ")
    res4,ws4 = OSSDP.solve(P,q,A,b,K,setMeanAdaptive);
    print("$(sr4.solverName): Cost: $(res4.cost), AvgTime/Iter: $(res4.avgIterTime) ")

    MeanData[kkk,:] = [res1.avgIterTime res2.avgIterTime res3.avgIterTime res4.avgIterTime]
    kkk+=1
    open(fn, "a") do f
      write(f,"ProblemType: $(pType) - ProblemName: $(problemName)\n")
      write(f, "$(sr1.solverName):\tCost: $(res1.cost),\t\tAvgTime/Iter: $(res1.avgIterTime)\n$(sr2.solverName):\tCost: $(res2.cost),\t\tAvgTime/Iter: $(res2.avgIterTime)\n$(sr3.solverName):\t\tCost: $(res3.cost),\t\tAvgTime/Iter: $(res3.avgIterTime)\n$(sr4.solverName):\t\tCost: $(res4.cost),\t\tAvgTime/Iter: $(res4.avgIterTime)\n\n")
    end

  end
  println(">>> ProblemType: $(pType) completed!")
end

# ---------------------------------------
# RUN SOME SDP LIB PROBLEMS
# ---------------------------------------
println("> RUN SDPLib Problems")

problemTypes = ["truss1";"hinf1"]
for problemName in problemTypes
  data = JLD.load("../sdplib/"*problemName*".jld")
  F = data["F"]
  c = data["c"]
  m = data["m"]
  n = data["n"]
  optVal = data["optVal"]

  P = zeros(m,m)
  q = c
  A = zeros(n^2,m)
  for iii = 1:m
    A[:,iii] = -vec(F[iii+1])
  end
  b = -vec(F[1])
  Kf = 0
  Kl = 0
  Kq = []
  Ks = [n^2]

  setUnNonAdaptive =   OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 0 ,scaleFunc=1,adaptive_rho=false)
  setMeanNonAdaptive = OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 10,scaleFunc=2,adaptive_rho=false)
  setUnAdaptive =   OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 0 ,scaleFunc=1,adaptive_rho=true)
  setMeanAdaptive = OSSDPSettings(max_iter=maxIter,checkTermination=1,scaling = 10,scaleFunc=2,adaptive_rho=true)

  # solve with QOCS
  K = Cone(Kf,Kl,Kq,Ks)

 # Solve with OSSDP
  println(problemName)
  res1,ws1 = OSSDP.solve(P,q,A,b,K,setUnNonAdaptive);
  print("QOCS - Un - noRho: Cost: $(res1.cost), AvgTime/Iter: $(res1.avgIterTime) ")
  res2,ws2 = OSSDP.solve(P,q,A,b,K,setMeanNonAdaptive);
  print("QOCS - Avg - noRho: Cost: $(res2.cost), AvgTime/Iter: $(res2.avgIterTime) ")
  res3,ws3 = OSSDP.solve(P,q,A,b,K,setUnAdaptive);
  print("QOCS - Un: Cost: $(res3.cost), AvgTime/Iter: $(res3.avgIterTime) ")
  res4,ws4 = OSSDP.solve(P,q,A,b,K,setMeanAdaptive);
  print("QOCS - Avg - Rho: Cost: $(res4.cost), AvgTime/Iter: $(res4.avgIterTime) ")

  MeanData[kkk,:] = [res1.avgIterTime res2.avgIterTime res3.avgIterTime res4.avgIterTime]
  kkk+=1
  open(fn, "a") do f
    write(f,"ProblemName: $(problemName)\n")
    write(f, "QOCS - Un - noRho:\tCost: $(res1.cost),\t\tAvgTime/Iter: $(res1.avgIterTime)\nQOCS - Avg - noRho:\tCost: $(res2.cost),\t\tAvgTime/Iter: $(res2.avgIterTime)\nQOCS - Un:\t\tCost: $(res3.cost),\t\tAvgTime/Iter: $(res3.avgIterTime)\nQOCS - Avg:\t\tCost: $(res4.cost),\t\tAvgTime/Iter: $(res4.avgIterTime)\n\n")
  end
end


# ---------------------------------------
# CALCULATE MEAN VALUE FOR ALL PROBLEMS
# ---------------------------------------
MeanData = MeanData[1:kkk-1,:]
open(fn, "a") do f
      write(f,"\nMean Data:\n")
      write(f, "QOCS - Un - noRho:\tMean: $(mean(MeanData[:,1]))\nQOCS - Avg - noRho:\tMean: $(mean(MeanData[:,2]))\nQOCS - Un:\tMean: $(mean(MeanData[:,3]))\nQOCS - Avg:\tMean: $(mean(MeanData[:,4]))\n")
  end



