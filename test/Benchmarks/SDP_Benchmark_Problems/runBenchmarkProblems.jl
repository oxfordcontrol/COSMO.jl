# Test routine that solve all SDP Benchmark problems and saves the results in .jld resultDataFiles
# uses the Compare Module

workspace()
include("../../../src/Solver.jl")
include("../solverComparison/Compare.jl")

using OSSDP, Base.Test, Compare, JLD
SAVE_ALWAYS = true
maxIter = 2000

# choose experiment name, otherwise save in folder based on timestamp
EXPERIMENT_NAME = "BadlyScaled"
# create a new folder for the results based on timestamp
timestamp = Dates.format(now(), "yyddmm_HH-MM")

isdefined(:EXPERIMENT_NAME) ? folderName=EXPERIMENT_NAME : folderName=timestamp
resPath = "../resultDataFiles/SDP_Benchmark_Problems/"*folderName
!ispath(resPath) && mkdir(resPath)

# find available problem types in SDP_Benchmark_Problems folder
# probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Julia/"
probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Scaled/"
existingFolders = readdir(probPath)
problemTypes = []
for f in filter(x -> !startswith(x, "."), readdir(probPath))
    f = split(f,".")[1]
    push!(problemTypes,String(f))
end
filter!(x->!in(x,["SmallestCircle"]),problemTypes)
println(">>> $(length(problemTypes)) Problem Type(s) detected!")

# run tests for each problem type
for pType in problemTypes
  println(">>> ProblemType: $(pType)")
  resFileName = resPath*"/"*pType*".jld"
  dirPath = probPath*pType*"/"

  # find all file names in problem type folder
  problems = []
  for f in filter(x -> endswith(x, ".jld"), readdir(dirPath))
      f = split(f,".")[1]
      push!(problems,String(f))
  end
  nn = length(problems)


  problemType = JLD.load("$(dirPath)"*"$(problems[1]).jld","problemType")
  # sr1 = SolverResult(nn, problemType,"COSMO - Un",timestamp,0,false,false)
  # sr2 = SolverResult(nn, problemType,"COSMO - Avg",timestamp,0,true,false)
  # sr3 = SolverResult(nn, problemType,"COSMO - Geom",timestamp,0,true,false)
  # sr4 = SolverResult(nn, problemType,"COSMO - Sym",timestamp,0,true,false)
  sr5 = SolverResult(nn, problemType,"COSMO - Un",timestamp,0,false,true)
  sr6 = SolverResult(nn, problemType,"COSMO - Avg",timestamp,0,true,true)
  sr7 = SolverResult(nn, problemType,"COSMO - Geo",timestamp,0,true,true)
  sr8 = SolverResult(nn, problemType,"COSMO - Sym",timestamp,0,true,true)
  # resData = [sr1;sr2;sr3;sr4;sr5;sr6;sr7;sr8]
  resData = [sr5;sr6;sr7;sr8]

  # loop over all problems in problem type folder
  for iii =1:1:length(problems)
    iii == length(problems) && (SAVE_ALWAYS = true)
    gc()
    problem = problems[iii]
    data = JLD.load("$(dirPath)"*"$(problem).jld")
    P = data["P"]
    q = data["q"]
    r = data["r"]
    A = data["A"]
    b = data["b"]
    m = data["m"]
    n = data["n"]
    Kf = data["Kf"]
    Kl = data["Kl"]
    Kq = data["Kq"]
    Ks = data["Ks"]

    # objTrue = data["objTrue"]
    problemName = data["problemName"]


    pDims = [size(A,1);size(A,2);nnz(A)]
    # update the true value for the QOCP solver
    # setUnNonAdaptive =   OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 0 ,scaleFunc=1,adaptive_rho=false)
    # setMeanNonAdaptive = OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 10,scaleFunc=2,adaptive_rho=false)
    # setGeoNonAdaptive =  OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 10,scaleFunc=3,adaptive_rho=false)
    # setSymNonAdaptive =  OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 10,scaleFunc=4,adaptive_rho=false)

    setUnAdaptive =   OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 0 ,scaleFunc=1,adaptive_rho=true)
    setMeanAdaptive = OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 10,scaleFunc=2,adaptive_rho=true)
    setGeoAdaptive =  OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 10,scaleFunc=3,adaptive_rho=true)
    setSymAdaptive =  OSSDPSettings(max_iter=maxIter,check_termination=1,scaling = 10,scaleFunc=4,adaptive_rho=true)

    # define cone membership
    K = Cone(Kf,Kl,Kq,Ks)


    # Solve with OSSDP
    # res1,nothing = OSSDP.solve(P,q,A,b,K,setUnNonAdaptive);
    # print("\n.")
    # res2,nothing = OSSDP.solve(P,q,A,b,K,setMeanNonAdaptive);
    # print(".")
    # res3,nothing = OSSDP.solve(P,q,A,b,K,setGeoNonAdaptive);
    # print(".")
    # res4,nothing = OSSDP.solve(P,q,A,b,K,setSymNonAdaptive);
    # print(".")
    res5,nothing = OSSDP.solve(P,q,A,b,K,setUnAdaptive);
    print(".")
    res6,nothing = OSSDP.solve(P,q,A,b,K,setMeanAdaptive);
    print(".")
    res7,nothing = OSSDP.solve(P,q,A,b,K,setGeoAdaptive);
    print(".")
    res8,nothing = OSSDP.solve(P,q,A,b,K,setSymAdaptive);
    print(".")

    # save results
    # resArray = [res1;res2;res3;res4;res5;res6;res7;res8]
    resArray = [res5;res6;res7;res8]
    updateResults!(resFileName,resData,resArray,pDims,problemName,r,SAVE_ALWAYS)
    printStatus(iii,nn,problemName,resData)
    # println("MOSEK sol: $(objTrue)")
  end
  println(">>> ProblemType: $(pType) completed!")
end

