# Unit tests for functions in src/algebra.jl file
using COSMO, Test, Random, LinearAlgebra, SparseArrays

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

@testset "Algebra" begin

for T in UnitTestFloats

  rng = Random.MersenneTwister(9191)


  @testset "Algebra (T = $(T))" begin

    # clip()
    n = 50
    x = 10 * rand(rng, T, n) .- 5
    x[10] = -8
    x[20] = 8
    x_c = zeros(T, n)
    x_c2 = zeros(T, n)
    @. x_c = COSMO.clip(x, -2, 2)
    @test (maximum(x_c) <= 2. && minimum(x_c) >= -2.)
    @. x_c2 = COSMO.clip(x, -2, 2, -100, 100)
    @test (x_c2[10] == -100 && x_c2[20] == 100)

    # scaled_norm
    E = Diagonal(T[1, 2, 3.])
    v = T[-1, 5, 4]
    @test COSMO.scaled_norm(E, v, 1) == 23
    @test COSMO.scaled_norm(E, v, 2) == norm(T[-1; 10; 12], 2)
    @test COSMO.scaled_norm(E, v, Inf) == 12
    @test_throws ArgumentError COSMO.scaled_norm(E, v, 3);

    # col_norms!()
    v = zeros(T, 4)
    A = T[1. 2 3 4; 2 -1 30 4.1]
    COSMO.col_norms!(v, A)
    @test v == T[2; 2; 30; 4.1]

    A = sparse(T[1. 2 3 4; 2 -1 30 4.1])
    B = sparse(T[1 3 -2 4.2; -100 -1 -2 -100])
    COSMO.col_norms!(v, A, reset = true)
    COSMO.col_norms!(v, B, reset = false)
    @test v == T[100; 3; 30; 100]

    # row_norms!()
    A = Matrix(A')
    B = Matrix(B')
    COSMO.row_norms!(v, A, reset = true)
    COSMO.row_norms!(v, B, reset = false)
    @test v == T[100; 3; 30; 100]
    A = sparse(A)
    B = sparse(B)
    COSMO.row_norms!(v, A, reset = true)
    COSMO.row_norms!(v, B, reset = false)
    @test v == T[100; 3; 30; 100]

    # lmul!()
    L = Diagonal(rand(rng, T, 5))
    Mo = sprand(rng, T, 5, 5, 0.8)
    M = deepcopy(Mo)
    COSMO.lmul!(L, M)
    @test norm(M -  L * Mo, Inf) <= 1e-10

    M = deepcopy(Mo)
    COSMO.lmul!(I, M)
    @test M == Mo

    x = rand(T, 5)
    y = copy(x)
    COSMO.lmul!(L, x)
    @test norm(x - L * y, Inf ) <= 1e-10

    # rmul!()
    R = Diagonal(rand(rng, T, 5))
    Mo = sprand(rng, T, 5, 5, 0.8)
    M = deepcopy(Mo)
    COSMO.rmul!(M, R)
    @test norm(M -  Mo * R, Inf) <= 1e-10

    M = deepcopy(Mo)
    COSMO.rmul!(M, I)
    @test M == Mo

    # lrmul!()
    L = Diagonal(rand(rng, T, 5))
    Mo = rand(rng, T, 5, 5)
    M = deepcopy(Mo)
    R = Diagonal(rand(rng, T, 5))
    COSMO.lrmul!(L, M, R)
    @test norm(M -  L * Mo * R, Inf) <= 1e-7

    M = deepcopy(Mo)
    COSMO.lrmul!(L, M, I)
    @test norm(M -  L * Mo, Inf) <= 1e-10

    M = deepcopy(Mo)
    COSMO.lrmul!(I, M, I)
    @test M == Mo

    # symmetrize!()
    A = rand(rng, T, 3, 3)
    B = deepcopy(A)
    COSMO.symmetrize_upper!(A)
    @test Symmetric(A) == (0.5 * (B + B'))
    @test_throws AssertionError COSMO.symmetrize_upper!(rand(T, 2, 3));

    # is_pos_def()
    A = Symmetric(rand(rng, T, 4, 4))
    Apos = A + 4 * diagm(0 => ones(T, 4))
    Aneg = A - 4 * diagm(0 => ones(T, 4))
    @test COSMO.is_pos_def(Apos, T(1e-6))
    @test !COSMO.is_pos_def(Aneg, T(1e-6))
    @test COSMO.is_neg_def(Aneg, T(1e-6))
    @test !COSMO.is_neg_def(Apos, T(1e-6))

  end
end
end
nothing
