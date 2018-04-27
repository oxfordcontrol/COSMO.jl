workspace()
include("../../src/Solver.jl")
include("../solverComparison/Compare.jl")
 include("./LatexExporter.jl")

# load data files
using Formatting, JLD, OSSDP, Compare
using LatexExport

timestamp = "182604_16-43"
dir = "../resultDataFiles/SDP_Benchmark_Problems/"*timestamp
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end



  metricsArr = Array{Compare.ProblemMetrics}(length(results))

# Step 1: loop over all problem types and the combined file and calculate important metrics
for iii=1:length(results)
  r = results[iii]
  data = JLD.load(dir*"/"*r*".jld")
  resData = data["resData"]

  contains(resData[1].problemType, "Combine") ? COMBINE_FLAG = true : COMBINE_FLAG = false
  # create data object to hold the metrics for this problem type
  pm = Compare.ProblemMetrics(resData[1].problemType,COMBINE_FLAG,resData[1].ind,length(resData))

  # loop over solver and compute metrics
  k = 1
  for s in resData
    solvedInd = find(x->x==:solved,s.status)
    # calculate mean
    meanIterAll = mean(s.iter)
    length(solvedInd) > 0 ? meanIterSolved = mean(s.iter[solvedInd]) : meanIterSolved = Inf

    numSolved = length(solvedInd)
    percSolved = numSolved/s.ind

    sm = Compare.SolverMetrics(s.solverName,s.adaptionON,s.scalingON,meanIterAll,meanIterSolved,numSolved,percSolved)
    pm.solverResults[k] = sm
    k+=1
  end
  metricsArr[iii] = pm
end

# permute metricsArr to have combined results at the end of array
cind = find(x-> x.combinedData,metricsArr)[1]
if cind == 1
  p = [collect(2:length(metricsArr));1]
elseif cind > 1 && cind < length(metricsArr)
  p = [collect(1:cind-1);collect(cind+1:length(metricsArr));2]
end
permute!(metricsArr,p)
# Step 2: Print results to screen
println("-"^80)
println("Postprocessing statistics:")
for pm in metricsArr
  println("-"^50)
  println(">>> Problem Type: $(pm.problemType)")
  println("- Combined Data: $(pm.combinedData)")
  println("- Number of Problems: $(pm.numProblems)")
  println("- Number of Solvers: $(pm.numSolvers)")

    println("Solver Name:\tAdaption ON:\tScaling ON:\tNum Solved:\t% Solved:\tMean Iter (all):\tMean Iter (solved):")
  for sm in pm.solverResults

    printfmt("{1:s}\t{2:b}\t\t{3:b}\t\t{4:d}\t\t{5:.2f}\t\t{6:.2f}\t\t\t{7:.2f}\n",sm.solverName,sm.adaptionON,sm.scalingON,sm.numSolved,sm.percSolved,sm.meanIterAll,sm.meanIterSolved)
    # println("Solver Name: $(sm.solverName[1:10])\tadaptionON: $(sm.adaptionON)\tscalingON: $(sm.scalingON)\tNum Solved: $(sm.numSolved)/$(pm.numProblems)\t% solved: $(sm.percSolved)\tMean Iter (all): $(sm.meanIterAll)\tMean Iter (solved): $(sm.meanIterSolved)")
    # println("Solver Name: $(sm.solverName)\t\tadaptionON: $(sm.adaptionON)\tscalingON: $(sm.scalingON)\tNum Solved: $(sm.numSolved)/$(pm.numProblems)\t% solved: $(sm.percSolved)\tMean Iter (all): $(sm.meanIterAll)\tMean Iter (solved): $(sm.meanIterSolved)")
  end
end
println("-"^80)

# Step 3: Print results to LaTeX table
resPath = dir*"/latex/"
!ispath(resPath) && mkdir(resPath)

 createLatexTable(metricsArr,resPath)



