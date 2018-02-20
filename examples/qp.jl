# Test script to test solver for a qp
workspace()
include("../src/Solver.jl")

using Base.Test
using OSSDP, OSSDPTypes

# Quadratic program example from OSQP Doc
# min 0.5 * x'Px +  q'x
# s.t. l <= x <= u
q = [1; 1]
P = sparse([4. 1; 1 2])
A = [1. 1; 1 0; 0 1]
l = [1.; 0; 0]
u = [1; 0.7; 0.7]

# create augmented matrices
Aa = [A eye(3) zeros(3,3);A zeros(3,3) -eye(3)]
ba = [u; l]
qa = [q;zeros(6,1)]
Pa = blkdiag(P,spzeros(3,3),spzeros(3,3))

# define cone (x1,x2, and 6 slack variables >= 0)
K = OSSDPTypes.Cone(2,6,[],[])
# define example problem
settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0)
res, = OSSDP.solve(Pa,qa,Aa,ba,K,settings)

@testset "QP Problem" begin
  @test isapprox(res.x[1:2],[0.3;0.7], atol=1e-3)
  @test isapprox(res.cost,1.88, atol=1e-3)
end

