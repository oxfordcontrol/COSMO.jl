# Test script to test solver for an lp

using COSMO, LinearAlgebra, SparseArrays, Test

# Linear program example
# min c'x
# s.t. Ax <= b
#      x >= 1,  x2 >= 5, x1+x3 >= 4
# with data matrices
c = [1; 2; 3; 4.]
A = Matrix(1.0I, 4, 4)
b = [10; 10; 10; 10]
n = 4
# -------------------
# create constraints A * x + b in set
# -------------------
# Ax <= b
c1 = COSMO.Constraint(-A, b, COSMO.Nonnegatives)
# x >= 1
c2 = COSMO.Constraint(Matrix(1.0I, n, n), -ones(n), COSMO.Nonnegatives)
# x2 >= 5
c3 = COSMO.Constraint(1, -5, COSMO.Nonnegatives, n, 2:2)
# x1 + x3 >= 4
c4 = COSMO.Constraint([1 0 1 0], -4, COSMO.Nonnegatives)

# -------------------
# define cost function
# -------------------
P = spzeros(4, 4)
q = c

# -------------------
# assemble solver model
# -------------------
settings = COSMO.Settings(max_iter=2500, verbose=true, eps_abs = 1e-4, eps_rel = 1e-5)
model = COSMO.Model()
assemble!(model, P, q, [c1; c2; c3; c4], settings = settings)
res = COSMO.optimize!(model);

@testset "Linear Problem" begin
  @test isapprox(res.x[1:4], [3; 5; 1; 1], atol=1e-2, norm = (x -> norm(x, Inf)))
  @test isapprox(res.obj_val, 20.0, atol=1e-2)
end
nothing
