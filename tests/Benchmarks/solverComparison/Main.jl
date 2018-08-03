workspace()
include("../../../src/Solver.jl")
include("./Compare.jl")

using OSQP, OSSDP, Compare,JLD

RUN_QP_LASSO_TEST = true
RUN_MESZAROS_TEST = true
RUN_SOCP_TEST = true
RUN_SDP_TEST = true
RUN_SVM_TEST = true
SAVE_ALWAYS = true

setOFF = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=1500,verbose=false,scaling = 0,check_termination = 1,eps_abs = 1e-3,eps_rel=1e-3,timelimit=60)
setON = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=1500,verbose=false,scaling = 10,check_termination = 1,eps_abs = 1e-3,eps_rel=1e-3,timelimit=60)


if RUN_MESZAROS_TEST
  timestamp = Dates.format(now(), "yyddmm_HH-MM")
  meszarosDataFile = "SC_" * timestamp * "meszarosQP.jld"
  dirPath = "../meszaros/JLD_Files/"
  problems = meszarosFilenames(dirPath)
  numProblems = length(problems)
  problemType = "QP-Meszaros"

  sr1 = SolverResult(numProblems, problemType,"OSSDP",timestamp,setOFF,false)
  sr2 = SolverResult(numProblems, problemType,"OSSDP",timestamp,setON,true)
  sr3 = SolverResult(numProblems, problemType,"OSQP",timestamp,0.,false)
  sr4 = SolverResult(numProblems, problemType,"OSQP",timestamp,0.,true)
  resData = [sr1;sr2;sr3;sr4]
  iii = 1
  for problem in problems
    # try
      #gc()
      # load problem data
      data = JLD.load("$(dirPath)"*"$(problem).jld")
      P1, q1,r1, A1, b1, K1 = loadMeszarosData(data,"QOCS")
      P2, q2,r2, A2, l2, u2 = loadMeszarosData(data,"OSQP")
      pDims = getMeszarosDim(data)

      # solve problem
      resOSSDP_unscaled,nothing = OSSDP.solve(P1,q1,A1,b1,K1,setOFF)
      print("\n.")
      resOSSDP_scaled,nothing = OSSDP.solve(P1,q1,A1,b1,K1,setON)
      print(".")

      m2u = OSQP.Model()
      OSQP.setup!(m2u; P=P2, q=q2, A=A2, l=l2, u=u2,scaling=0,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
      tic()
      resOSQP_unscaled = OSQP.solve!(m2u)
      sr3.runTime[sr3.ind+1] = toq()
      print(".")

      m2s = OSQP.Model()
      OSQP.setup!(m2s; P=P2, q=q2, A=A2, l=l2, u=u2,scaling=10,max_iter=1500,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
      tic()
      resOSQP_scaled = OSQP.solve!(m2s)
      sr4.runTime[sr4.ind+1] = toq()
      print(".\n")

      # update and save result
      resArray = [resOSSDP_unscaled;resOSSDP_scaled;resOSQP_unscaled;resOSQP_scaled]
      updateResults!(meszarosDataFile,resData,resArray,pDims,problem,r1,SAVE_ALWAYS)
      printStatus(iii,numProblems,problem,resData)
      iii +=1
    # catch
    #   println("An error happened with problem $(problem).")
    #   continue
    # end
  end

end