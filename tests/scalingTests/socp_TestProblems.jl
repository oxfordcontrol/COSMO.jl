# Test routine to compare scaling for a number of SOCP Lasso problems (partially badly scaled)
workspace()
include("../../src/Solver.jl")

using OSSDP, Base.Test, JLD

rng = MersenneTwister(12345)

nn = 25

resIter = zeros(nn,6)
resCost = zeros(nn,6)



for iii =1:1:nn
  # generate problem data
  n = 10
  m = 100*n
  F = rand(m,n)

  vtrue = sprand(n,1,0.1 )
  noise = 0.1*rand(m,1)
  b = F*vtrue + noise
  μMax = norm(F'*b,Inf)
  μ = 0.1*μMax

  # create augmented matrices
  Aa = [eye(n) eye(n) zeros(n,1) -eye(n) zeros(n,n) zeros(n,1) zeros(n,m+1);
      -eye(n) eye(n)  zeros(n,1) zeros(n,n) eye(n)  zeros(n,1) zeros(n,m+1);
      zeros(1,n) zeros(1,n) -1   zeros(1,n) zeros(1,n) 1 zeros(1,m+1);
      zeros(1,n) zeros(1,n) 1   zeros(1,n) zeros(1,n) 0 1 zeros(1,m);
      zeros(m,n) -2*F zeros(m,1) zeros(m,n) zeros(m,n) zeros(m,1) zeros(m,1) eye(m,m)]

  ba = [zeros(n,1);zeros(n,1);1;1;-2*b]
  ca = [μ*ones(n,1);zeros(n,1);0.5;zeros(2*n+1+1+m)]
  Pa = zeros(length(ca),length(ca))

  # define cone (t, z and w free variables, 2*n slack variables for inequality constraints, m+2 soc variables)
  K = Cone(2*n+1,2*n,[2+m],[])

  # Solve with OSSDP
   settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled100,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=1000.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled1000,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled100,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
  settings = OSSDPSettings(rho=1000.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled1000,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)


  resCost[iii,:] = [resOSSDP_unscaled.cost resOSSDP_scaled.cost resOSSDP_unscaled100.cost resOSSDP_scaled100.cost resOSSDP_unscaled1000.cost resOSSDP_scaled1000.cost]
  resIter[iii,:] = [resOSSDP_unscaled.iter resOSSDP_scaled.iter resOSSDP_unscaled100.iter resOSSDP_scaled100.iter resOSSDP_unscaled1000.iter resOSSDP_scaled1000.iter]
  println("$(iii)/$(nn) completed!")
end


timestamp = Dates.format(now(), "yyddmm_HH-MM")
filename = timestamp * "_socpScalingTest.jld"
JLD.save(filename, "resCost", resCost, "resIter",resIter)
println(">>> Test Data successfully saved in $(filename).\n")
