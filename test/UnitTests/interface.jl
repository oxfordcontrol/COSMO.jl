using COSMO, Test, LinearAlgebra, Random
rng = Random.MersenneTwister(41)


# set!() -  to directly set the problem data without building constraints first

P = [4. 1;1 2]
q = [1; 1.]
A = [1. 1;1 0; 0 1]
l = [1.;0;0]
u = [1.;0.7;0.7]
Aa = [A; -A]
b = [u; -l]
convex_sets = [COSMO.Nonnegatives(3); COSMO.Nonnegatives(3)]

model = COSMO.Model()
COSMO.set!(model, P, q, Aa, b, convex_sets )
res = COSMO.optimize!(model)
@test isapprox(res.obj_val, 1.88, atol = 1e-3)

qf = [1; 1; 1]
Pf = zeros(1, 1)
Af = [1 2; 1 2]
bf = [1; 2]
model = COSMO.Model()
@test_throws DimensionMismatch COSMO.set!(model, P, qf, Aa, b, convex_sets)
@test_throws DimensionMismatch COSMO.set!(model, Pf, q, Aa, b, convex_sets)
@test_throws DimensionMismatch COSMO.set!(model, P, q, Af, b, convex_sets)
@test_throws DimensionMismatch COSMO.set!(model, P, q, Aa, bf, convex_sets)

