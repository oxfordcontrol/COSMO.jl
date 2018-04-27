module LatexExport

using Compare
export createLatexTable


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

end #MODULE