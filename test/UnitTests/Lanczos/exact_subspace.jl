using COSMO, Test, LinearAlgebra, SparseArrays, Random

# generate problem data
rng = Random.MersenneTwister(12345)
n = 100 

X0 = Symmetric(randn(rng, n, n))
@testset "Lanczos Projection - Exact Subspace" begin
  for offset in [-10, 10]
      X = X0 + offset*I;
      # @show sum(eigen(X).values .> 0)
      cone = COSMO.PsdCone(n^2)
      x1 = X[:]
      COSMO.project_exact!(x1, cone)
      x2 = X[:]
      COSMO.project!(x2, cone)
      @test isapprox(x1, x2, atol=1e-2, norm=(x -> norm(x,Inf)))
      # Test numerical stability of the subspace
      for i = 1:30
        x2 = X[:]
        COSMO.project!(x2, cone)
      end
      @test isapprox(x1, x2, atol=1e-2, norm=(x -> norm(x,Inf)))
  end
end

nothing