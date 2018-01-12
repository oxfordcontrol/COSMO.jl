# Test script to solve a easy QP for debugging purposses
workspace()
include("../Projections.jl")
include("../Solver.jl")

using Base.Test
using OSSDP
# using PyPlot


#-------------------------
# LP Case:
# ------------------------

# Problem Data
P = zeros(2,2)
Q = [1.0 0;0 -1]
q = vec(Q)
Av = [0.0 0;0 1]
A = [0.0 0 0 1]
b = 5.0
# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.0,max_iter=30,verbose=true)

# solve SDP problem
res,dbg = solveSDP(P,q,A,b,settings);
res

#-------------------------
# QP Case: 0.5 x^2 - 3x
# ------------------------

# QP Problem Data
P = 1.0
q = -3.0
A = 0.0
b = 0.0
# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true)

# solve SDP problem
res,dbg = solveSDP(P,q,A,b,settings)

@testset "QP: 0.5 x^2 - 3x results" begin
  @test norm(res.x[1]- 3.0) < 1e-3
  @test norm(res.cost +4.5) < 1e-3
  @test norm(res.x-res.s,Inf)<1e-3
  @test res.x[1] > 0.0
end


#-------------------------
# QP Case: 0.5 x^2 + 5x
# ------------------------

# QP Problem Data
P = 1.0
q = 5.0
A = 0.0
b = 0.0
# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.0,max_iter=2500,verbose=true)

# solve SDP problem
res,dbg = solveSDP(P,q,A,b,settings)

@testset "QP: 0.5 x^2 + 5x results" begin
  @test norm(res.x[1]- 0) < 1e-3
  @test norm(res.cost ) < 1e-3
  @test norm(res.x-res.s,Inf)<1e-3
  @test res.x[1] >= 0.0
end

# # compute smallest eigenvalue of result X
# F = eigfact(X)
# XλMin = minimum(F[:values])
# F = eigfact(reshape(res.s,3,3))
# SλMin = minimum(F[:values])
# epsEV = -1e-15
# @testset "Compare results" begin
#   @test norm(X-[1.05593 0.369185 0.868298; 0.369185 0.129079 0.303585; 0.868299 0.303585 0.714011],Inf) < 1e-3
#   @test norm(res.cost-13.902255839911007) < 1e-3
#   @test norm(res.x-res.s,Inf)<1e-3
#   @test XλMin > epsEV
#   @test SλMin > epsEV
# end


