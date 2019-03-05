using COSMO, Test, LinearAlgebra, Random, SparseArrays
rng = Random.MersenneTwister(41)


@testset "Interface" begin
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


# assemble! function test different combinations of P an q

#P not sparse, q dense
P = rand(rng, 2, 2)
q = rand(rng, 2)
A = rand(rng, 1, 2)
b = 1.
constraint = COSMO.Constraint(A, b, COSMO.Nonnegatives)
model = COSMO.Model();
COSMO.assemble!(model, P, q, constraint)
@test model.p.P == P
@test model.p.q == q
@test model.p.A == -A
@test model.p.b == [b]

# sparse input data
COSMO.empty_model!(model)
COSMO.assemble!(model, sparse(P), sparse(q), constraint)
@test model.p.P == P
@test model.p.q == q

A = rand(rng, 5, 1)
b = rand(rng, 5)
constraint = COSMO.Constraint(A, b, COSMO.Nonnegatives)

# P number, q vector
COSMO.empty_model!(model)
COSMO.assemble!(model, 1., [1.], constraint)
@test model.p.P == hcat([1.])
@test model.p.q == [1.]

# P vector, q number
COSMO.empty_model!(model)
COSMO.assemble!(model, [1.], 1., constraint)
@test model.p.P == hcat([1.])
@test model.p.q == [1.]

# P number, q number
COSMO.empty_model!(model)
COSMO.assemble!(model, 1., 1., constraint, settings = COSMO.Settings(verbose = true))
@test model.p.P == hcat([1.])
@test model.p.q == [1.]

# P vector, q matrix
COSMO.empty_model!(model)
COSMO.assemble!(model, [1.], hcat([1.]), constraint)
@test model.p.P == hcat([1.])
@test model.p.q == [1.]

# Test case where dimension of an A in one constraint is inconsistent with P
A = 1.
b = 0.
P = sparse(1.0I, 2, 2)
q = rand(rng, 2)
constraint = COSMO.Constraint(A, b, COSMO.Nonnegatives)
model = COSMO.Model();
@test_throws DimensionMismatch COSMO.assemble!(model, P, q, [constraint])

end