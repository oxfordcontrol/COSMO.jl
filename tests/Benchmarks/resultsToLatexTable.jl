# File that saves test results into a latex table

using PyPlot, JLD

data = JLD.load("./180402_12-47_sdpLibResults_scaling.jld")
resX = data["resX"]
resCost = data["resCost"]
resIter = data["resIter"]
resRT = data["resRT"]
problems = data["problems"]
resStatus = data["status"]
resPrim = data["resPrim"]
resDual = data["resDual"]
set = data["settings"]
resIter[:,1] = [6;7;7;8;7;11;17]
resIter[:,2] = [140;140;1180;3320;8560;20000;20000]



n = length(problems)
# 1. plot errors in objective value for scs vs mosek and ossdp vs mosek for different problems
errors = zeros(n,3)
errors[:,1] = abs(resCost[:,1]-resCost[:,2])
errors[:,2] = abs(resCost[:,1]-resCost[:,3])
errors[:,3] = abs(resCost[:,1]-resCost[:,4])


caption1 = "Comparison of MOSEK, SCS and OSSDP"

open("./latex/table.tex", "w") do file
  # create the header of the table
  write(file, "\\begin{longtable}{c c c c c c}\n\\toprule\n")
  write(file, "  &  & MOSEK & SCS & OSSDP & OSSDP (scaled) \\\\ [0.5ex]\n\\hline\n")
  # add data
  iii = 1
  for p in problems
    write(file, "$(p) & obj & $(round(resCost[iii,1],4))& $(round(resCost[iii,2],4)) & $(round(resCost[iii,3],4)) & $(round(resCost[iii,4],4)) \\\\ \n & error & & $(round(errors[iii,1],6)) & $(round(errors[iii,2],6)) & $(round(errors[iii,3],6)) \\\\ \n  & rPrim & &  & $(round(resPrim[iii,3],6)) & $(round(resPrim[iii,4],6)) \\\\ \n& rDual & & & $(round(resDual[iii,3],6)) & $(round(resDual[iii,4],6)) \\\\ \n& iter & $(resIter[iii,1]) &$(resIter[iii,2]) & $(resIter[iii,3]) & $(resIter[iii,4]) \\\\ \n & status & $(resStatus[iii,1]) &$(resStatus[iii,2]) & $(resStatus[iii,3]) & $(resStatus[iii,4]) \\\\ \n \\hline\n")
    iii +=1
  end
  # write footer
  write(file, "\\bottomrule\n\\caption{$(caption1)}\n\\label{table:1}\n\\end{longtable}")
  # print out the solver parameters as well
  write(file, "\n \n \\begin{table}[h]\n\\begin{tabular}{l l |l l |l l}\n")
  write(file, " \$\\sigma\$ & $(set.sigma) & \$\\rho\$ & $(set.rho) & \$\\alpha\$ & $(set.alpha) \\\\ \n\\hline\n \$\\epsilon_{abs}\$ & $(set.eps_abs) & \$\\epsilon_{rel}\$ & $(set.eps_rel) & scale & $(set.scaling) \\\\ \n\\hline\n checkTerm & $(set.checkTermination) &  & \\\\\n\\hline\n")
  write(file, "\n\\end{tabular}\n\\centering\n\\caption{Parameters used in OSSDP Solver}\n\\end{table}")


end





