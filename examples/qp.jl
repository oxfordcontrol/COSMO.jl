# Test script to test solver for a qp

using Test
using QOCS, SparseArrays, LinearAlgebra

# Quadratic program example from OSQP Doc
# min 0.5 * x'Px +  q'x
# s.t. l <= Ax <= u
q = [1; 1]
P = sparse([4. 1; 1 2])
A = [1. 1; 1 0; 0 1]
l = [1.; 0; 0]
u = [1; 0.7; 0.7]

# create augmented matrices
Aa = [A;-A]
ba = [u; -l]

# define cone (x1,x2, and 6 slack variables >= 0)
K = QOCS.Cone(0,6,[],[])
# define example problem
settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,check_termination=1,adaptive_rho = true, scaling = 10,eps_abs = 1e-4,eps_rel = 1e-4)
res, ws= QOCS.solve(P,q[:],Aa,ba[:],K,settings)

@testset "QP Problem" begin
  @test norm(res.x[1:2] - [0.3;0.7],Inf) < 1e-3
  @test abs(res.cost-1.88) < 1e-3
end
nothing
