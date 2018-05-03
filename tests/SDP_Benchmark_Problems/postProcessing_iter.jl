# File to determine in how many cases a scaled solver was better than a unscaled solver (only for cases with adaptive rho)

workspace()
include("../../src/Solver.jl")
include("../solverComparison/Compare.jl")
 include("./LatexExporter.jl")

# load data files
using Formatting, JLD, OSSDP, Compare
using LatexExport

folderName = "AfterBugFix"
dir = "../resultDataFiles/SDP_Benchmark_Problems/"*folderName
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end
# filter!(x->!in(x,["Combined"]),results)

# permute to have combined results at the end of array
cind = find(x-> x == "Combined",results)
if length(cind) > 0
  cind = cind[1]
  if cind == 1
    p = [collect(2:length(results));1]
  elseif cind > 1 && cind < length(results)
    p = [collect(1:cind-1);collect(cind+1:length(results));2]
  end
  permute!(results,p)
end

resComp = zeros(Int64,length(results),4)
numProblemsArr = zeros(length(results))
println("-"^80)
println("Postprocessing statistics:")
# Step 1: loop over all problem types and the combined file and calculate important metrics
for iii=1:length(results)
  r = results[iii]
  data = JLD.load(dir*"/"*r*".jld")
  resData = data["resData"]
  numProblems = resData[1].ind
  numProblemsArr[iii] = numProblems
  unscaledIter = resData[5].iter
  avg = resData[6].iter
  geo = resData[7].iter
  sym = resData[8].iter
  bestScaled = minimum([avg geo sym],2)

  resComp[iii,1] = length(filter(x->x, avg .<= unscaledIter))
  resComp[iii,2] = length(filter(x->x, geo .<= unscaledIter))
  resComp[iii,3] = length(filter(x->x, sym .<= unscaledIter))
  resComp[iii,4] = length(filter(x->x, bestScaled .<= unscaledIter))

  println("-"^50)
  println(">>> Problem Type: $(r), Number of Problems: $(numProblems)")
  println("Avg <= Unscaled: $(resComp[iii,1]) ($(round(resComp[iii,1]/numProblems*100,2)) %)")
  println("Geo <= Unscaled: $(resComp[iii,2]) ($(round(resComp[iii,2]/numProblems*100,2)) %)")
  println("Sym <= Unscaled: $(resComp[iii,3]) ($(round(resComp[iii,3]/numProblems*100,2)) %)")
  println("Any Scaled <= Unscaled: $(resComp[iii,4]) ($(round(resComp[iii,4]/numProblems*100,2)) %)")

end
println("-"^80)


# Step 3: Print results to LaTeX table
resPath = dir*"/latex/"
!ispath(resPath) && mkdir(resPath)
 # createIterLatexTable(resComp,numProblemsArr,results,resPath)



