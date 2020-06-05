# Unit Test for constraint type

using COSMO, Test, SparseArrays, Random

# This test is precision agnostic and will be run with TestFloat precision
if @isdefined UnitTestFloat
    TestFloat = UnitTestFloat
else
    TestFloat = Float64
end

rng = Random.MersenneTwister(1872381)

A_int = 4
b_int = 2

A_Float = TestFloat(4.)
b_Float = TestFloat(2.)

A_UInt = UInt64(4)
b_UInt = UInt64(2)

A_sparse = sprand(rng, TestFloat, 10, 2, 0.4)
b_sparse = sparse(rand(rng, TestFloat, 10))

A_rvec = TestFloat[1 2 3 4]
b_rvec = TestFloat(1)

A_colvec = TestFloat[1; 2; 3; 4]
b_colvec = TestFloat[4; 3; 2; 1]

A_mat = rand(rng, TestFloat, 10, 10)
b_mat = sparse(rand(rng, TestFloat, 10, 1))



@testset "Constraints" begin

  @testset "Constructors" begin
    @test typeof(COSMO.Constraint(A_int, b_int, COSMO.ZeroSet)) <: COSMO.Constraint{Float64} #integers are always converted to Float64
    @test typeof(COSMO.Constraint(A_Float, b_Float, COSMO.ZeroSet)) <: COSMO.Constraint{TestFloat}
    @test typeof(COSMO.Constraint(A_UInt, b_UInt, COSMO.ZeroSet)) <: COSMO.Constraint{Float64}
    @test typeof(COSMO.Constraint(A_sparse, b_sparse, COSMO.ZeroSet)) <: COSMO.Constraint{TestFloat}
    @test typeof(COSMO.Constraint(A_rvec, b_rvec, COSMO.ZeroSet)) <: COSMO.Constraint{TestFloat}
    @test typeof(COSMO.Constraint(A_colvec, b_colvec, COSMO.ZeroSet)) <: COSMO.Constraint{TestFloat}
    @test typeof(COSMO.Constraint(A_mat, b_mat, COSMO.ZeroSet)) <: COSMO.Constraint{TestFloat}
  end

  @testset "Indices" begin
    A = rand(rng, TestFloat, 3, 3)
    b = rand(rng, TestFloat, 3)
    ind = 3:5
    dim = 10
    cs = COSMO.Constraint(A, b, COSMO.ZeroSet, dim, ind)

    @test cs.A[:, ind] == A
   end

  @testset "Merge Constraints" begin
    A1 = rand(rng, TestFloat, 10, 10)
    b1 = rand(rng, TestFloat, 10)
    A2 = rand(rng, TestFloat, 10, 10)
    b2 = rand(rng, TestFloat, 10)
    Am = [A1; A2]
    bm = [b1; b2]
    c1 = COSMO.Constraint(A1, b1, COSMO.ZeroSet)
    c2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)
    cArr = [c1; c2]
    cm = COSMO.Constraint(Am, bm, COSMO.ZeroSet)
    COSMO.merge_constraints!(cArr)

    @test cArr[1].A == cm.A
    @test cArr[1].b == cm.b
    @test cArr[1].convex_set.dim == cm.convex_set.dim
    @test typeof(cArr[1].convex_set) == typeof(cm.convex_set)


    c1 = COSMO.Constraint(A1, b1, COSMO.Nonnegatives)
    c2 = COSMO.Constraint(A2, b2, COSMO.Nonnegatives)
    cArr = [c1; c2]
    cm = COSMO.Constraint(Am, bm, COSMO.Nonnegatives)
    COSMO.merge_constraints!(cArr)

    @test cArr[1].A == cm.A
    @test cArr[1].b == cm.b
    @test cArr[1].convex_set.dim == cm.convex_set.dim
    @test typeof(cArr[1].convex_set) == typeof(cm.convex_set)
  end





end
nothing
