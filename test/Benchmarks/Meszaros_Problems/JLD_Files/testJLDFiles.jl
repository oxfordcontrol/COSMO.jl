workspace()

using OSQP, JLD, Base.Test, MAT

fileNames = []
for f in filter(x -> endswith(x, ".jld"), readdir("."))
    f = split(f,".")[1]
    push!(fileNames,String(f))
end
readmeInfo = JLD.load("../MAT_FILES/objVals.jld")
problemData = readmeInfo["problemData"]
sortedInd = sort!(collect(1:1:length(fileNames)), by=i->problemData[i,4])
fileNames = fileNames[sortedInd]
fileNames = fileNames[1:135]
iii = 1
@testset "Test JLD files" begin
  for file in fileNames
      # load problem data
      data = JLD.load("./"*"$(file).jld")
      P = data["P"]
      A = data["A"]
      q = data["q"]
      u = data["u"]
      l = data["l"]
      r = data["r"]

      m1 = OSQP.Model()
      OSQP.setup!(m1; P=P, q=q[:], A=A, l=l[:], u=u[:],scaling=10,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
      res = OSQP.solve!(m1)
      solvedCostJLD = res.info.obj_val + r


      data2 = MAT.matread("../bart_meszaros_data/"*"$(file).mat")
      P2 = data2["P"]
      A2 = data2["A"]
      q2 = data2["q"]
      u2 = data2["u"]
      l2 = data2["l"]
      r2 = data2["r"]

      m2 = OSQP.Model()
      OSQP.setup!(m2; P=P2, q=q2[:], A=A2, l=l2[:], u=u2[:],scaling=10,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
      res2 = OSQP.solve!(m2)
      solvedCostMAT = res2.info.obj_val + r2




       # check result
      objTrue = data["objVal"]

      @test (abs(solvedCostJLD-solvedCostMAT) < 1e-2)

      println("$(iii)/$(length(fileNames)): $(abs(solvedCostJLD-solvedCostMAT))")
      iii+=1
  end

end