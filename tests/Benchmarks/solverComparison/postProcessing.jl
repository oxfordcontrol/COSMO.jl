workspace()
include("../../../src/Solver.jl")

include("./Compare.jl")

using PyPlot, JLD, Compare, OSSDPTypes, Base.Test
cc =  ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b"]


data = JLD.load("../resultDataFiles/SC_181803_18-04MaxEigSDP.jld")
resData = data["resData"]
nn = resData[1].ind
eTol = 1e-2


# nameValDict = JLD.load("./../meszaros/JLD_Files/objVals.jld")["nameValDict"]
# optVal = Dict()
# for iii = 1:resData[1].ind
#   name = lowercase(resData[1].problemName[iii])
#   optVal[resData[1].problemName[iii]] = nameValDict[name]
# end

# function findSolvedInd(resData,k)
#   solInd = []
#   iii = 1
#   other = collect(1:1:length(resData))
#   iseven(k) ? filter!(x->x!=k && x!=(k-1),other) : filter!(x->x!=k && x!=(k+1),other)

#   for p in resData[1].problemName[1:resData[1].ind]
#     ck = resData[k].objVal[iii]
#     cArr = map(x->resData[x].objVal[iii],other)
#     objTrue = optVal[p]
#     if abs(cArr[1]-ck) < 1e-2 || abs(cArr[2]-ck) < 1e-2  || abs( ck-objTrue) < 1e-2
#       push!(solInd,iii)
#     end
#     iii +=1
#   end
#   return solInd
# end



# #----------------------
# # Print results summary
# #----------------------
# iii = 1
# myind =0
# println("Number of problems: $(nn)")
# for s in resData
#   indSolved = find(x->(x == :solved || x == :Solved || x == :Solved_inaccurate),s.status)
#   indMaxIter = find(x->(x == :UserLimit || x == :Max_iter_reached),s.status)
#   indTimelimit = find(x->(x == :TimeLimit),s.status)
#   indCloseToSol = find(x->abs(s.objVal[x]-optVal[s.problemName[x]]) < eTol,collect(1:s.ind))
#   indClose = findSolvedInd(resData,iii)
#   s.scalingON ? (si = "scaled") : (si = "unscaled")
#   println("Solver: $(s.solverName) - $(si), Solved: $(length(indSolved)), CloseToTrueSol: $(length(indCloseToSol)), CloseToSol: $(length(indClose)), MaxIter: $(length(indMaxIter)), TimeLimit: $(length(indTimelimit))")
#   iii +=1
# end



# #----------------------
# # Did scaling improve convergence?
# #----------------------
# println("----------- Did scaling improve convergence? -----------")
# for iii in [1 3]
#     indCloseOFF = findSolvedInd(resData,iii)
#     indCloseON = findSolvedInd(resData,iii+1)
#     commonInd = union(indCloseON,indCloseOFF)
#     # find out if iterations were reduced by scaling the data
#     OFFwins = commonInd[find(x->resData[iii].iter[x] < resData[iii+1].iter[x],commonInd)]
#     ONwins = commonInd[find(x->resData[iii].iter[x] > resData[iii+1].iter[x],commonInd)]
#     tie = commonInd[find(x->resData[iii].iter[x] == resData[iii+1].iter[x],commonInd)]
#     # find out how many iterations were saved or added on average
#     OFFmeanSavedIter = mean(resData[iii+1].iter[OFFwins]-resData[iii].iter[OFFwins])
#     ONmeanSavedIter = mean(resData[iii].iter[ONwins]-resData[iii+1].iter[ONwins])
#     # print out results
#     println("Solver: $(resData[iii].solverName): Better without scaling: $(length(OFFwins)), Better with scaling: $(length(ONwins)), Tie: $(length(tie)), Avg saved iter when off wins:$(OFFmeanSavedIter), Avg saved iter when on wins: $(ONmeanSavedIter)")

# end




#----------------------
# create iteration plot for all problems (solved or not)
#----------------------
PyPlot.figure(1,facecolor="white",figsize=(12,5))
iii = 1
for s in resData
  !s.scalingON ? ls = "--" : ls="-"
  s.adaptionON ? ad = "- adapt" : ad =""
  s.scalingON ? (si = "scaled") : (si = "unscaled")
  contains(s.solverName,"QOCS") && !s.adaptionON && (c = cc[1])
  contains(s.solverName,"QOCS") && s.adaptionON && (c = cc[2])
  contains(s.solverName,"OSQP") && !s.adaptionON && (c = cc[3])
  contains(s.solverName,"OSQP") && s.adaptionON && (c = cc[4])
  PyPlot.plot(1:1:s.ind,s.iter[1:s.ind], ls,color=c,label=s.solverName*" - $(si)$(ad)")
  iii+=1
end
PyPlot.grid(true)
PyPlot.xticks(1:1:nn)
PyPlot.ylabel("Iterations to Convergence",fontsize=15)
PyPlot.legend(ncol = 4,bbox_to_anchor=(0., 1.02, 1., .102),loc=3)
PyPlot.xticks(1:1:nn, resData[1].problemName, size="small",rotation="vertical")
#PyPlot.savefig("./iter.eps")


# #----------------------
# # create iteration plot only for problems solved by all solvers (configurations)
# #----------------------
# commonInd = findSolvedInd(resData,1)
# for iii=2:length(resData)
#     indCloseToSol = findSolvedInd(resData,iii)
#     commonInd = intersect(commonInd,indCloseToSol)
# end

# PyPlot.figure(2,facecolor="white",figsize=(12,5))
# iii = 1
# Ncom = length(commonInd)
# for s in resData
#   !s.scalingON ? ls = "--" : ls="-"
#   s.scalingON ? (si = "scaled") : (si = "unscaled")
#   contains(s.solverName,"OSSDP") ? c = cc[1] : c = cc[2]
#   PyPlot.plot(1:1:Ncom,s.iter[commonInd], ls,color=c,label=s.solverName*" - $(si)")
#   iii+=1
# end
# PyPlot.grid(true)
# PyPlot.xticks(1:1:Ncom)
# PyPlot.xlabel("Problem",fontsize=15)
# PyPlot.ylabel("Iterations to Convergence",fontsize=15)
# PyPlot.legend(ncol = 4,bbox_to_anchor=(0., 1.02, 1., .102),loc=3)
# PyPlot.title(resData[1].problemType*" - solved problems")
# PyPlot.xticks(1:1:Ncom, resData[1].problemName[commonInd], size="small",rotation="vertical")
# nothing
