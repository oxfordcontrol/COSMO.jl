# Test script to solve a easy QP for debugging purposses
workspace()
include("../Projections.jl")
include("../Solver.jl")

using Base.Test
using OSSDP
# using PyPlot


#-------------------------
# QP Case
# ------------------------

# QP Problem Data
P = 1
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]
# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true)

# solve SDP problem
res,dbg = solveSDP(P,q,A,b,settings)
X = reshape(res.x,3,3)


# compute smallest eigenvalue of result X
F = eigfact(X)
XλMin = minimum(F[:values])
F = eigfact(reshape(res.s,3,3))
SλMin = minimum(F[:values])
epsEV = -1e-15
@testset "Compare results" begin
  @test norm(X-[1.05593 0.369185 0.868298; 0.369185 0.129079 0.303585; 0.868299 0.303585 0.714011],Inf) < 1e-3
  @test norm(res.cost-13.902255839911007) < 1e-3
  @test norm(res.x-res.s,Inf)<1e-3
  @test XλMin > epsEV
  @test SλMin > epsEV
end


