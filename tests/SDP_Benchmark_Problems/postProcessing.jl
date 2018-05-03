# workspace()
# include("../../src/Solver.jl")
# include("../solverComparison/Compare.jl")

# using PyPlot, JLD, Compare, OSSDPTypes, Base.Test
cc =  ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b"]

folderName = "182604_10-19"
dir = "../resultDataFiles/SDP_Benchmark_Problems/"*folderName
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end

# MAIN LOOP
for iii=1:length(results)
  r = results[iii]
  data = JLD.load(dir*"/"*r*".jld")
  resData = data["resData"]
  # FIXME: ONLY TEMPORARY CHANGE!!!!!!!!!!
  resData = resData[1:3]
  problemName = resData[1].problemName[1]
  nn = resData[1].ind
  eTol = 1e-2


  PyPlot.figure(iii,facecolor="white",figsize=(12,5))
  kkk = 1
  for s in resData
    !s.scalingON ? ls = "--" : ls="-"
    s.adaptionON ? ad = "- adapt" : ad =""
    # s.scalingON ? (si = "scaled") : (si = "unscaled")
    # contains(s.solverName,"QOCS") && !s.adaptionON && (c = cc[1])
    # contains(s.solverName,"QOCS") && s.adaptionON && (c = cc[2])
    if kkk == 1
      c = cc[1]
      si = "unscaled"
    elseif kkk == 2
      c = cc[2]
      si = "Ruiz (Avg)"
    else
      c = cc[3]
      si = "Ruiz (Sym)"
    end
    PyPlot.plot(1:1:s.ind,s.iter[1:s.ind], ls,color=c,label=s.solverName*" - $(si)$(ad)")
    kkk+=1
  end
  PyPlot.grid(true)
  PyPlot.xticks(1:1:nn)
  PyPlot.title(resData[1].problemType, loc="left")
  PyPlot.xlabel("Problem number",fontsize=15)
  PyPlot.ylabel("Iterations to Convergence",fontsize=15)
  PyPlot.legend(ncol = 2,bbox_to_anchor=(0., 1.05, 1., .102))
  # PyPlot.xticks(1:1:nn, resData[1].problemName, size="small",rotation="vertical")
  PyPlot.savefig(dir*"/$(timestamp)-$(problemName).eps")
end


