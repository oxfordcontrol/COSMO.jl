# Unit Test for by split vector type
using COSMO, Test, Random, LinearAlgebra

rng = Random.MersenneTwister(13131)

convex_set = COSMO.PsdConeTriangle(3)
composite_set = COSMO.CompositeConvexSet([COSMO.ZeroSet(3), COSMO.Nonnegatives(1)])
x1 = rand(3)
x2 = rand(4)

s1 = COSMO.SplitVector(x1, convex_set)
s2 = COSMO.SplitVector(x2, composite_set)
@test size(s1) == (3, )
@test length(s2) == 4
@test getindex(s1, 2) == x1[2]
setindex!(s1, 5, 2)
@test getindex(s1, 2) == 5.
@test firstindex(s2) == 1
@test lastindex(s2) == 4
@test Base.showarg(IOBuffer(), s2, true ) == nothing