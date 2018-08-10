  # Test routine to compare scaling for a number of SDP problems (partially badly scaled)
workspace()
include("../../../src/Solver.jl")

using OSSDP, Base.Test

rng = Random.MersenneTwister(12345)

nn = 20

resIter = zeros(20,2)
resCost = zeros(20,2)


for iii =1:1:nn
  # generate problem data
  m = 1
  r = 8
  n = r^2
  A = Symmetric(100*sprand(r,r,0.5))
  A = vec(A)'
  # create strictly feasible primal point
  X = Symmetric(10*sprand(r,r,0.5))
  X = X + (-minimum(eigs(X)[1])+1)*speye(r)
  b = A*vec(X)
  b = [b]

  # create strictly feasibile dual point
  y = rand(m,1)
  S = Symmetric(10*sprand(r,r,0.5))
  S = S + (-minimum(eigs(S)[1])+1)*speye(r)
  c = vec(S) + A'*y

  P = zeros(n,n)

  K = Cone(0,0,[],[n])

  # Solve with OSSDP
  settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,check_termination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_unscaled,nothing = OSSDP.solve(P,c,A,b,K,settings)
  settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,check_termination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
  resOSSDP_scaled,nothing = OSSDP.solve(P,c,A,b,K,settings)

  resCost[iii,:] = [resOSSDP_unscaled.cost resOSSDP_scaled.cost]
  resIter[iii,:] = [resOSSDP_unscaled.iter resOSSDP_scaled.iter]
  println("$(iii)/$(nn) completed!")
end


# @testset "Lasso Scaling" begin
#     @test abs(resOSSDP_unscaled.cost-resOSQP1_unscaled.info.obj_val) < 1e-2
#     @test norm(resOSSDP_unscaled.x[1:n]-resOSQP1_unscaled.x[1:n],Inf) < 1e-2
#     @test norm(resOSSDP_unscaled.x[1:n]-resOSQP2_unscaled.x[1:n],Inf) < 1e-2
#     @test norm(resOSQP1_unscaled.x[1:n]-resOSQP2_unscaled.x[1:n],Inf) < 1e-2
#   end