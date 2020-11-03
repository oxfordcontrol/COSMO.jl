# SOCP Lasso Testproblem

using COSMO, Test, LinearAlgebra, SparseArrays, Random


# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

@testset "SOCP - Lasso" begin
  for T in UnitTestFloats
      if T != BigFloat
      # generate problem data
      rng = Random.MersenneTwister(12345)
      n = 8
      m = 50 * n
      F = rand(rng, T, m, n)

      vtrue = sprand(rng, T, n, 1, 0.1)
      noise = T(0.1) * rand(rng, T, m, 1)
      b = F * vtrue + noise
      μMax = norm(F' * b, Inf)
      μ = T(0.1) * μMax


      # define lasso problem as SOCP

      A1 = -sparse([1 zeros(T, 1, 2 * n + 1) 1 zeros(T, 1, m);
        -1 zeros(T, 1, 2 * n) 1 zeros(T, 1, m + 1);
        zeros(T, m, 1) -2 * F zeros(T, m, n + 2) diagm( 0 => ones(T, m))])

      A2 = -sparse([zeros(T, n, 1) diagm(0 => ones(T, n)) -diagm(0 => ones(T, n)) zeros(T, n, m + 2);
        zeros(T, n, 1) -diagm(0 => ones(T, n)) -diagm(0 => ones(T, n)) zeros(T, n, m + 2)])
      A3 = -sparse([zeros(T, 1, 2 * n + 1) -one(T) zeros(T, 1, m + 1);
       zeros(T, 1, 2 * n + 2) -one(T) zeros(T, 1, m);
       zeros(T, m, 2 * n + 3) -diagm( 0 => ones(T, m))])
      b1 = T[1; 1; -2 * b]
      b2 = zeros(T, 2 * n)
      b3 = zeros(T, m + 2)

      q = [one(T); zeros(T, n); μ * ones(T, n); zeros(T, m + 2)]
      P = spzeros(T, length(q), length(q))

      cs1 = COSMO.Constraint(A1 ,b1, COSMO.ZeroSet)
      cs2 = COSMO.Constraint(A2 ,b2, COSMO.Nonnegatives)
      cs3 = COSMO.Constraint(A3 ,b3, COSMO.SecondOrderCone)


      # Solve with OSSDP
      model = COSMO.Model{T}()
      assemble!(model, P, q, [cs1; cs2; cs3])
      res = COSMO.optimize!(model);

      @testset "SOCP - Lasso (T = $(T))" begin
      @test res.status == :Solved
      @test isapprox(res.obj_val, 0.4422849814458825,atol=1e-2)
      end
    end
  end
end
nothing
