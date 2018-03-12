# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
workspace()
include("../../src/Solver.jl")

using OSQP, OSSDP, Base.Test, JLD

rng = MersenneTwister(12345)

nn = 25

resIter = zeros(nn,8)
resCost = zeros(nn,8)


for iii =1:1:nn
  # generate problem data
  n = 10
  m = 100*n
  A = sprandn(rng,m,n,0.5)
  vtrue = 1/n*sprandn(rng,n,0.5)
  noise = 1/4*randn(rng,m,1)
  b = A*vtrue + noise
  λ = 0.2*norm(A'*b,Inf)

  # define lasso problem as QP
  Aa = [-A zeros(m,n) eye(m,m) zeros(m,2*n);
         eye(n,n) -eye(n,n) zeros(n,m) eye(n,n) zeros(n,n);
         eye(n,n) eye(n,n) zeros(n,m) zeros(n,n) -eye(n,n) ]

  ba = [b;zeros(2*n)]
  P = 2*diagm([zeros(2*n);ones(m);zeros(2*n)]) # times two to cancel the 1/2 in the cost function
  q = [zeros(n);λ*ones(n);zeros(m+2*n)]

  # define cone membership
  K = OSSDPTypes.Cone(2*n+m,2*n,[],[])

  # Solve with OSSDP
  settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled,nothing = OSSDP.solve(P,q,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled,nothing = OSSDP.solve(P,q,Aa,ba,K,settings)

  settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled100,nothing = OSSDP.solve(P,q,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=1000.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled1000,nothing = OSSDP.solve(P,q,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled100,nothing = OSSDP.solve(P,q,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=1000.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled1000,nothing = OSSDP.solve(P,q,Aa,ba,K,settings)


  # modify problem for OSQP (using inequality constraints)
  Aa1 = sparse(Aa[1:2*n+m,1:m+2*n])
  l = full([b;-Inf*ones(n);zeros(n)][:])
  u = full([b;zeros(n);Inf*ones(n)][:])
  q1 = q[1:2*n+m]
  P1 = sparse(P[1:2*n+m,1:m+2*n])
  m1 = OSQP.Model()
  m1s = OSQP.Model()
  OSQP.setup!(m1; P=P1, q=q1, A=Aa1, l=l, u=u,scaling=0,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
  OSQP.setup!(m1s; P=P1, q=q1, A=Aa1, l=l, u=u,scaling=10,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
  resOSQP1_unscaled = OSQP.solve!(m1)
  resOSQP1_scaled = OSQP.solve!(m1s)


  # # modify problem for OSQP (using equality constraints)
  # Aa2 = sparse([Aa;zeros(n,2*n+m) eye(n,n) zeros(n,n);zeros(n,2*n+m) zeros(n,n) ones(n,n) ])
  # l2 = [l;zeros(2*n)];
  # u2 = [u;Inf*ones(2*n)];
  # m2 = OSQP.Model()
  # m2s = OSQP.Model()

  # OSQP.setup!(m2; P=sparse(P), q=q, A=Aa2, l=l2, u=u2,scaling=0,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
  # OSQP.setup!(m2s; P=sparse(P), q=q, A=Aa2, l=l2, u=u2,scaling=10,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
  # resOSQP2_unscaled = OSQP.solve!(m2)
  # resOSQP2_scaled = OSQP.solve!(m2s)

  resCost[iii,:] = [resOSSDP_unscaled.cost resOSSDP_scaled.cost resOSSDP_unscaled100.cost resOSSDP_scaled100.cost resOSSDP_unscaled1000.cost resOSSDP_scaled1000.cost resOSQP1_unscaled.info.obj_val resOSQP1_scaled.info.obj_val]# resOSQP2_unscaled.info.obj_val resOSQP2_scaled.info.obj_val]
  resIter[iii,:] = [resOSSDP_unscaled.iter resOSSDP_scaled.iter resOSSDP_unscaled100.iter resOSSDP_scaled100.iter resOSSDP_unscaled1000.iter resOSSDP_scaled1000.iter resOSQP1_unscaled.info.iter resOSQP1_scaled.info.iter]# resOSQP2_unscaled.info.iter resOSQP2_scaled.info.iter]
  println("$(iii)/$(nn) completed!")
end

timestamp = Dates.format(now(), "yyddmm_HH-MM")
filename = timestamp * "_qpScalingTest.jld"
JLD.save(filename, "resCost", resCost, "resIter",resIter)
println(">>> Test Data successfully saved in $(filename).\n")

