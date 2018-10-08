# Test routine to compare performance for a number of SOCP Lasso problems
workspace()
include("../../../src/Solver.jl")
include("../solverComparison/Compare.jl")

using OSQP, OSSDP, Base.Test, Compare, JuMP, Mosek

rng = Random.MersenneTwister(12345)

nn = 25
timestamp = Dates.format(now(), "yyddmm_HH-MM")
dataFile = "./SC_" * timestamp * "LassoSOCP.jld"
problemType = "SOCP-Lasso"

sr1 = SolverResult(nn, problemType,"COSMO",timestamp,0,false)
sr2 = SolverResult(nn, problemType,"COSMO",timestamp,0,true)
resData = [sr1;sr2]
ws = 0
for iii =1:1:nn
  # generate problem data
  n = 8
  m = 50*n
  F = rand(m,n)

  vtrue = sprand(n,1,0.1 )
  noise = 0.1*rand(m,1)
  b = F*vtrue + noise
  μMax = norm(F'*b,Inf)
  μ = 0.1*μMax

  pDims = [m;n;m*n]

  # define lasso problem as SOCP

  # create augmented matrices
  Aa = sparse([1 zeros(1,2*n+1) 1 zeros(1,m);
        -1 zeros(1,2*n) 1 zeros(1,m+1);
        zeros(m,1) -2*F zeros(m,n+2) Matrix(1.0I,m,m);
        zeros(n,1) eye(n) -eye(n) zeros(n,m+2);
        zeros(n,1) -eye(n) -eye(n) zeros(n,m+2);
       zeros(1,2*n+1) -1 zeros(1,m+1);
       zeros(1,2*n+2) -1 zeros(1,m);
       zeros(m,2n+3) -eye(m)])

  ba = [1;1;-2*b;zeros(2*n+m+2)]
  q = [1;zeros(n);μ*ones(n,1);zeros(m+2,1)]
  P = spzeros(length(q),length(q))

  # define cone membership
  K = Cone(2+m,2*n,[m+2],[])

setOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,check_termination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-3,adaptive_rho=true)
setON = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,check_termination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-3,adaptive_rho=true)

  # Solve with OSSDP
  resOSSDP_unscaled,nothing = OSSDP.solve(P,q,Aa,ba,K,setOFF);
  print("\n.")
  resOSSDP_scaled,nothing = OSSDP.solve(P,q,Aa,ba,K,setON);
  print(".")

  model = Model(solver=MosekSolver())
  @variable(model, x[1:n])
  @variable(model, t[1:n])
  @variable(model, y)
  @variable(model, w[1:m+1])
  @objective(model, Min, y + μ*ones(n)'*t)
  @constraint(model, -t .<= x)
  @constraint(model, x .<= t)
  @constraint(model, w[1] == 1-y)
  @constraint(model, w[2:m+1] .== 2*(F*x-b))
  @constraint(model, norm(w) <= 1+y)
  status = JuMP.solve(model)
  println("Objective value: ", getobjectivevalue(model))
  println("x = ", getvalue(x))

  resArray = [resOSSDP_unscaled;resOSSDP_scaled]
  updateResults!(dataFile,resData,resArray,pDims,"SOCP-$(iii)",0.,true)
  printStatus(iii,nn,"SOCP-$(iii)",resData)
end


