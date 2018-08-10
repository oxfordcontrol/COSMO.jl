workspace()
include("../../../src/Solver.jl")
include("../meszaros/ConvertProblem.jl")

include("./Compare.jl")

using PyPlot, JLD, Compare, OSSDPTypes, Base.Test
cc =  ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b"]


data1 = JLD.load("./SC_181603_15-22LassoQP.jld")
data2 = JLD.load("./SC_181503_10-04LassoQP.jld")
resData1 = data1["resData"]
resData2 = data2["resData"]
resData2[1].solverName = "DiffRho"
resData2[2].solverName = "DiffRho"
resData = [resData1;resData2]
nn = resData[1].ind
eTol = 1e-2

PyPlot.figure(1,facecolor="white",figsize=(12,5))
iii = 1
for s in resData
  !s.scalingON ? ls = "--" : ls="-"
  s.scalingON ? (si = "scaled") : (si = "unscaled")
  contains(s.solverName,"OSSDP") ? c = cc[1] : c = cc[2]
  PyPlot.plot(1:1:s.ind,s.iter[1:s.ind], ls,color=c,label=s.solverName*" - $(si)")
  iii+=1
end
PyPlot.grid(true)
PyPlot.xticks(1:1:nn)
PyPlot.xlabel("Problem",fontsize=15)
PyPlot.ylabel("Iterations to Convergence",fontsize=15)
PyPlot.legend(ncol = 4,bbox_to_anchor=(0., 1.02, 1., .102),loc=3)
PyPlot.title(resData[1].problemType*" - all problems")
PyPlot.xticks(1:1:nn, resData[1].problemName, size="small",rotation="vertical")