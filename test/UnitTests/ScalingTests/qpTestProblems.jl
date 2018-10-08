# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
workspace()
include("../../../src/Solver.jl")
include("../../Benchmarks/solverComparison/Compare.jl")

using OSQP, OSSDP, Base.Test, Compare, JuMP, Mosek

rng = Random.MersenneTwister(12345)
#rng = Random.MersenneTwister(211121)

nn = 20
timestamp = Dates.format(now(), "yyddmm_HH-MM")
dataFile = "./SC_" * timestamp * "LassoQP.jld"
problemType = "QP-Lasso"

sr1 = SolverResult(nn, problemType,"COSMO",timestamp,0,false)
sr2 = SolverResult(nn, problemType,"COSMO",timestamp,0,true)
sr3 = SolverResult(nn, problemType,"OSQP",timestamp,0,false)
sr4 = SolverResult(nn, problemType,"OSQP",timestamp,0,true)
resData = [sr1;sr2;sr3;sr4]
ws = 0
for iii =1:1:nn
  # generate problem data
  n = 8
  m = 50*n
  A = sprandn(rng,m,n,0.5)
  vtrue = 1/n*sprandn(rng,n,0.5)
  noise = 1/4*randn(rng,m,1)
  b = A*vtrue + noise
  λ = 0.2*norm(A'*b,Inf)

  pDims = [m;n;nnz(A)]

  # define lasso problem as QP
  Aa = [-A zeros(m,n) Matrix(1.0I,m,m);
         Matrix(1.0I,n,n) -Matrix(1.0I,n,n) zeros(n,m);
         -Matrix(1.0I,n,n) -Matrix(1.0I,n,n) zeros(n,m)]

  ba = [-b;zeros(2*n)]
  P = 2*Matrix(Diagonal([zeros(2*n);ones(m)])) # times two to cancel the 1/2 in the cost function
  q = [zeros(n);λ*ones(n);zeros(m)]

  # define cone membership
  K = OSSDPTypes.Cone(m,2*n,[],[])

  setOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,check_termination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-3,adaptive_rho=true)
  setON = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,check_termination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-3,adaptive_rho=true)

  # Solve with OSSDP
  resOSSDP_unscaled,nothing = OSSDP.solve(P,q,Aa,ba,K,setOFF);
        print("\n.")

  resOSSDP_scaled,nothing = OSSDP.solve(P,q,Aa,ba,K,setON);
      print(".")

  # modify problem for OSQP (using inequality constraints)
  l = full([-b;-Inf*ones(n);zeros(n)][:])
  u = full([-b;zeros(n);Inf*ones(n)][:])
  Aa2 = [-A zeros(m,n) Matrix(1.0I,m,m);
       Matrix(1.0I,n,n) -Matrix(1.0I,n,n) zeros(n,m);
       Matrix(1.0I,n,n) Matrix(1.0I,n,n) zeros(n,m)]

  m1 = OSQP.Model()
  OSQP.setup!(m1; P=sparse(P), q=q, A=sparse(Aa2), l=l, u=u,scaling=0,check_termination=1,verbose=false,adaptive_rho = true,eps_abs = 1e-3,eps_rel=1e-3)
  resOSQP1_unscaled = OSQP.solve!(m1)
  m2 = OSQP.Model()
  OSQP.setup!(m2; P=sparse(P), q=q, A=sparse(Aa2), l=l, u=u,scaling=10,check_termination=1,verbose=false,adaptive_rho = true,eps_abs = 1e-3,eps_rel=1e-3)
  resOSQP1_scaled = OSQP.solve!(m2)

  resArray = [resOSSDP_unscaled;resOSSDP_scaled;resOSQP1_unscaled;resOSQP1_scaled]
  updateResults!(dataFile,resData,resArray,pDims,"QP-$(iii)",0.,true)
  printStatus(iii,nn,"QP-$(iii)",resData)
end


