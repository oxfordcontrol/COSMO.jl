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

# Define the constraint l <= Ax <= u with the help of a Nonnegatives set
Aa = [-A;A]
ba = [u; -l]
constraint1 = QOCS.Constraint(Aa,ba,QOCS.Nonnegatives())

# define example problem
settings = QOCS.Settings(verbose=true,eps_abs = 1e-4,eps_rel = 1e-4)


model = QOCS.Model()
assemble!(model,P,q,constraint1)
res = QOCS.optimize!(model,settings);

# solve again by defining the constraints with the help of a box (disable infeasibility checks)
constraint1 = QOCS.Constraint(A,zeros(3),QOCS.Box(l,u))
settings = QOCS.Settings(check_infeasibility = 2500, verbose=true,eps_abs = 1e-4,eps_rel = 1e-4)

model = QOCS.Model()
assemble!(model,P,q,constraint1)
res_box = QOCS.optimize!(model,settings);


@testset "QP Problem" begin
  @test norm(res.x[1:2] - [0.3;0.7],Inf) < 1e-3
  @test norm(res_box.x[1:2] - [0.3;0.7],Inf) < 1e-3
  @test abs(res.objVal-1.88) < 1e-3
  @test abs(res_box.objVal-1.88) < 1e-3
end
nothing
