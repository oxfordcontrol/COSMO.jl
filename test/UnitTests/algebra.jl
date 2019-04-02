# Unit tests for functions in src/algebra.jl file
using COSMO, Test, Random, LinearAlgebra, SparseArrays

rng = Random.MersenneTwister(9191)


@testset "Algebra" begin

  # clip()
  n = 50
  x = 10 * rand(rng, n) .- 5
  x[10] = -8
  x[20] = 8
  x_c = zeros(n)
  x_c2 = zeros(n)
  @. x_c = COSMO.clip(x, -2, 2)
  @test (maximum(x_c) <= 2. && minimum(x_c) >= -2.)
  @. x_c2 = COSMO.clip(x, -2, 2, -100, 100)
  @test (x_c2[10] == -100 && x_c2[20] == 100)

  # scaled_norm
  E = Diagonal([1, 2, 3.])
  v = [-1, 5, 4]
  @test COSMO.scaled_norm(E, v, 1) == 23
  @test COSMO.scaled_norm(E, v, 2) == norm([-1; 10; 12], 2)
  @test COSMO.scaled_norm(E, v, Inf) == 12
  @test_throws ArgumentError COSMO.scaled_norm(E, v, 3);

  # col_norms!()
  v = zeros(4)
  A = [1. 2 3 4; 2 -1 30 4.1]
  COSMO.col_norms!(v, A)
  @test v == [2; 2; 30; 4.1]

  A = sparse([1. 2 3 4; 2 -1 30 4.1])
  B = sparse([1 3 -2 4.2; -100 -1 -2 -100])
  COSMO.col_norms!(v, A, reset = true)
  COSMO.col_norms!(v, B, reset = false)
  @test v == [100; 3; 30; 100]

  # row_norms!()
  A = Matrix(A')
  B = Matrix(B')
  COSMO.row_norms!(v, A, reset = true)
  COSMO.row_norms!(v, B, reset = false)
  @test v == [100; 3; 30; 100]
  A = sparse(A)
  B = sparse(B)
  COSMO.row_norms!(v, A, reset = true)
  COSMO.row_norms!(v, B, reset = false)
  @test v == [100; 3; 30; 100]

  # lmul!()
  L = Diagonal(rand(rng, 5))
  Mo = sprand(rng, 5, 5, 0.8)
  M = deepcopy(Mo)
  COSMO.lmul!(L, M)
  @test norm(M -  L * Mo, Inf) <= 1e-10

  M = deepcopy(Mo)
  COSMO.lmul!(I, M)
  @test M == Mo

  # rmul!()
  R = Diagonal(rand(rng, 5))
  Mo = sprand(rng, 5, 5, 0.8)
  M = deepcopy(Mo)
  COSMO.rmul!(M, R)
  @test norm(M -  Mo * R, Inf) <= 1e-10

  M = deepcopy(Mo)
  COSMO.rmul!(M, I)
  @test M == Mo

  # lrmul!()
  L = Diagonal(rand(rng, 5))
  Mo = rand(rng, 5, 5)
  M = deepcopy(Mo)
  R = Diagonal(rand(rng, 5))
  COSMO.lrmul!(L, M, R)
  @test norm(M -  L * Mo * R, Inf) <= 1e-10

  M = deepcopy(Mo)
  COSMO.lrmul!(L, M, I)
  @test norm(M -  L * Mo, Inf) <= 1e-10

  M = deepcopy(Mo)
  COSMO.lrmul!(I, M, I)
  @test M == Mo

  # symmetrize!()
  A = rand(rng, 3, 3)
  B = deepcopy(A)
  COSMO.symmetrize_upper!(A)
  @test Symmetric(A) == (0.5 * (B + B'))
  @test_throws AssertionError COSMO.symmetrize_upper!(rand(2, 3));

  # is_pos_sem_def()
  A = Symmetric(randn(rng, 4, 4))
  Apos = A + 4 * Matrix(1.0I, 4, 4)
  Aneg = A - 4 * Matrix(1.0I, 4, 4)
  @test COSMO.is_pos_sem_def(Apos, 1e-6)
  @test !COSMO.is_pos_sem_def(Aneg, 1e-6)
  @test COSMO.is_neg_sem_def(Aneg, 1e-6)
  @test !COSMO.is_neg_sem_def(Apos, 1e-6)

end
nothing