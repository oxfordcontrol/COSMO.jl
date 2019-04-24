# Test script to test solver for a qp

using Test
using COSMO, SparseArrays, LinearAlgebra

# Quadratic program example from OSQP Doc
# min 0.5 * x'Px +  q'x
# s.t. l <= Ax <= u
q = [1; 1.]
P = sparse([4. 1; 1 2])
A = [1. 1; 1 0; 0 1]
l = [1.; 0; 0]
u = [1; 0.7; 0.7]

# Define the constraint l <= Ax <= u with the help of a Nonnegatives set
Aa = [-A; A]
ba = [u; -l]
constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives)

# define example problem
settings = COSMO.Settings(verbose=true, eps_abs = 1e-4, eps_rel = 1e-4)


model = COSMO.Model()
assemble!(model, P, q, constraint1, settings = settings)
res = COSMO.optimize!(model);

# solve again by defining the constraints with the help of a box
constraint1 = COSMO.Constraint(A, zeros(3), COSMO.Box(l, u))

model = COSMO.Model()
assemble!(model, P, q, constraint1, settings = settings)
res_box = COSMO.optimize!(model);


@testset "QP Problem" begin
  @test norm(res.x[1:2] - [0.3; 0.7], Inf) < 1e-3
  @test norm(res_box.x[1:2] - [0.3; 0.7], Inf) < 1e-3
  @test abs(res.obj_val - 1.88) < 1e-3
  @test abs(res_box.obj_val - 1.88) < 1e-3
end
nothing
