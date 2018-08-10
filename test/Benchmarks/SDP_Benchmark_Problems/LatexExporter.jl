module LatexExport

using Compare
export createLatexTable, createIterLatexTable


function createLatexTable(metricsArr::Array{Compare.ProblemMetrics},dir::String)


  caption1 = "Performance of different solver configurations over the set of SDP Benchmark problems."

  open(dir*"/table.tex", "w") do file
    # create the header of the table
    write(file, "\\begin{longtable}{l c c c c c c}\n\\toprule\n")
    # add data
    for pm in metricsArr
      write(file, "\\multicolumn{7}{l}{Problem Type: \\textbf{$(pm.problemType)}} \\\\ \n \\multicolumn{7}{l}{Number of Problems: $(pm.numProblems)} \\\\ \n \\multicolumn{7}{l}{Number of Solvers: $(pm.numSolvers)} \\\\ \n  \\hline\n")
      write(file, " Solver Name &  \$\\rho\$-Adap & Scaling & \\# solved & \\% solved & Mean Iter(all) & Mean Iter(solved) \\\\ [0.5ex]\n\\hline\n")
      for sm in pm.solverResults
        write(file, "$(sm.solverName) & $(sm.adaptionON) & $(sm.scalingON) & $(sm.numSolved) & $(round(sm.percSolved,2)) & $(round(sm.meanIterAll,2)) & $(round(sm.meanIterSolved,2))  \\\\ \n")
      end
      write(file,"\\bottomrule \\\\ \n")
    end
    # write footer
    write(file, "\\bottomrule \n\\caption{$(caption1)}\n\\label{table:1}\n\\end{longtable}")

  end
  return nothing
end

  function createIterLatexTable(resComp,numProblemsArr,results,dir)
  caption1 = "Number of problems where a solver configuration with scaling performed better than the unscaled configuration (only cases with rho-adaption)."
  open(dir*"/iterTable.tex", "w") do file
    write(file, "\\begin{longtable}{l c c c }\n\\toprule\n")
    for iii=1:length(numProblemsArr)
      numProblems = numProblemsArr[iii]
      write(file, "\\multicolumn{4}{l}{Problem Type: \\textbf{$(results[iii])}} \\\\ \n \\multicolumn{4}{l}{Number of Problems: $(numProblemsArr[iii])} \\\\ \n  \\hline\n")
      write(file, " Mean \$ \\leq \$ Unscaled &  Geom \$ \\leq \$ Unscaled & Sym \$ \\leq \$ Unscaled &  Any Scaled \$ \\leq \$  Unscaled \\\\ [0.5ex]\n")
      write(file, "$(resComp[iii,1]) ($(round(resComp[iii,1]/numProblems*100,2)) \\%) & $(resComp[iii,2]) ($(round(resComp[iii,2]/numProblems*100,2)) \\%) & $(resComp[iii,3]) ($(round(resComp[iii,3]/numProblems*100,2)) \\%)  & $(resComp[iii,4]) ($(round(resComp[iii,4]/numProblems*100,2)) \\%)  \\\\ \n")
      write(file,"\\bottomrule \\\\ \n")
    end
        write(file, "\\bottomrule \n\\caption{$(caption1)}\n\\label{table:2}\n\\end{longtable}")

  end
  return nothing

  end
end #MODULE