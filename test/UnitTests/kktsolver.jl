using COSMO, Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(2401)

@testset "refactoring checks" begin

  m = 20
  n = 30

  P  = sparse(generate_pos_def_matrix(n,rng))
  A  = sprandn(m,n,0.2)
  sigma = rand(n)
  rho1  = rand(m)
  rho2  = rand(m)
  b = randn(m+n)


  F = COSMO.QdldlKKTSolver(P,A,sigma,rho1)
  J = [P+Diagonal(sigma) A'; A -Diagonal(1 ./ rho1)]
  x = copy(b)

  #test a single solve
  COSMO.solve!(F,x)
  @test norm(x - J\b) <= 1e-10

  #test a rho update and solve
  J = [P+Diagonal(sigma) A'; A -Diagonal(1 ./ rho2)]
  x = copy(b)
  COSMO.update_rho!(F,rho2)
  COSMO.solve!(F,x)
  @test norm(x - J\b) <= 1e-10

end
