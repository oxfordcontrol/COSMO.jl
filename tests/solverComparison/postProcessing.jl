workspace()
include("../../src/Solver.jl")
include("../meszaros/ConvertProblem.jl")

include("./Compare.jl")

using PyPlot, JLD, Compare, OSSDPTypes, Base.Test
cc =  ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b"]

data = JLD.load("./SC_181103_13-43meszarosQP.jld")
resData = data["resData"]
nn = resData[1].ind
eTol = 1e-2

nameValDict = JLD.load("./../meszaros/MAT_FILES/objVals.jld")["nameValDict"]
optVal = Dict()
for iii = 1:resData[1].ind
  name = lowercase(resData[1].problemName[iii])
  optVal[resData[1].problemName[iii]] = nameValDict[name]
end
#----------------------
# Print results summary
#----------------------
iii = 1
myind =0
println("Number of problems: $(nn)")
for s in resData
  indSolved = find(x->(x == :solved || x == :Solved || x == :Solved_inaccurate),s.status)
  indMaxIter = find(x->(x == :UserLimit || x == :Max_iter_reached),s.status)
  indTimelimit = find(x->(x == :TimeLimit),s.status)
  indCloseToSol = find(x->abs(s.objVal[x]-optVal[s.problemName[x]]) < eTol,collect(1:s.ind))
  println("Solver: $(s.solverName), Solved: $(length(indSolved)), CloseToSol $(length(indCloseToSol)), MaxIter: $(length(indMaxIter)), TimeLimit: $(length(indTimelimit))")
  iii +=1
end

#----------------------
# Did scaling improve convergence?
#----------------------
println("----------- Did scaling improve convergence? -----------")
for iii in [1 3]
    indCloseOFF = find(x->abs(resData[iii].objVal[x]-optVal[resData[iii].problemName[x]]) < eTol,collect(1:resData[iii].ind))
    indCloseON = find(x->abs(resData[iii+1].objVal[x]-optVal[resData[iii+1].problemName[x]]) < eTol,collect(1:resData[iii+1].ind))
    commonInd = union(indCloseON,indCloseOFF)
    # find out if iterations were reduced by scaling the data
    OFFwins = commonInd[find(x->resData[iii].iter[x] < resData[iii+1].iter[x],commonInd)]
    ONwins = commonInd[find(x->resData[iii].iter[x] > resData[iii+1].iter[x],commonInd)]
    tie = commonInd[find(x->resData[iii].iter[x] == resData[iii+1].iter[x],commonInd)]
    # find out how many iterations were saved or added on average
    OFFmeanSavedIter = mean(resData[iii+1].iter[OFFwins]-resData[iii].iter[OFFwins])
    ONmeanSavedIter = mean(resData[iii].iter[ONwins]-resData[iii+1].iter[ONwins])
    # ADD CODE HERE
    # print out results
    println("Solver: $(split(resData[iii].solverName,"-")[1]): Better without scaling: $(length(OFFwins)), Better with scaling: $(length(ONwins)), Tie: $(length(tie)), Avg saved iter when off wins:$(OFFmeanSavedIter), Avg saved iter when on wins: $(ONmeanSavedIter)")

end

#----------------------
# create iteration plot for all problems (solved or not)
#----------------------
PyPlot.figure(1,facecolor="white",figsize=(12,5))
iii = 1
for s in resData
  contains(s.solverName,"unscaled") ? ls = "--" : ls="-"
  contains(s.solverName,"OSSDP") ? c = cc[1] : c = cc[2]
  PyPlot.plot(1:1:s.ind,s.iter[1:s.ind], ls,color=c,label=s.solverName)
  iii+=1
end
PyPlot.grid(true)
PyPlot.xticks(1:1:nn)
PyPlot.xlabel("Problem",fontsize=15)
PyPlot.ylabel("Iterations to Convergence",fontsize=15)
PyPlot.legend(ncol = 4,bbox_to_anchor=(0., 1.02, 1., .102),loc=3)
PyPlot.title(resData[1].problemType*" - all problems")
PyPlot.xticks(1:1:nn, resData[1].problemName, size="small",rotation="vertical")
PyPlot.savefig("./iter.eps")


#----------------------
# create iteration plot only for problems solved by all solvers (configurations)
#----------------------
commonInd = find(x->abs(resData[1].objVal[x]-optVal[resData[1].problemName[x]]) < eTol,collect(1:resData[1].ind))
for iii=2:length(resData)
    s = resData[iii]
    indCloseToSol = find(x->abs((s.objVal[x]-optVal[s.problemName[x]])) < eTol,collect(1:s.ind))
    commonInd = intersect(commonInd,indCloseToSol)
end
PyPlot.figure(2,facecolor="white",figsize=(12,5))
iii = 1
Ncom = length(commonInd)
for s in resData
  contains(s.solverName,"unscaled") ? ls = "--" : ls="-"
  contains(s.solverName,"OSSDP") ? c = cc[1] : c = cc[2]
  PyPlot.plot(1:1:Ncom,s.iter[commonInd], ls,color=c,label=s.solverName)
  iii+=1
end
PyPlot.grid(true)
PyPlot.xticks(1:1:Ncom)
PyPlot.xlabel("Problem",fontsize=15)
PyPlot.ylabel("Iterations to Convergence",fontsize=15)
PyPlot.legend(ncol = 4,bbox_to_anchor=(0., 1.02, 1., .102),loc=3)
PyPlot.title(resData[1].problemType*" - solved problems")
PyPlot.xticks(1:1:Ncom, resData[1].problemName[commonInd], size="small",rotation="vertical")
nothing
