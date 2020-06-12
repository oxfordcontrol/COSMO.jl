# Dual infeasible test problems
# Choose A, P and q in a way such that the problem becomes unbounded below
# here the last element of x appears in the cost function with a negative sign and is unconstrained


# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

for T in UnitTestFloats
  if precision(T) >= precision(Float64)

    rng = Random.MersenneTwister(555)
    nn = 1

    @testset "Dual infeasible QP problems - Testset 1" begin
      for iii = 1:nn

      # choose size of problem
      n = rand(rng,5:50)
      m = 2*n
      A = sprand(rng, T, m,n,0.7)*T(50)
      A[:,end] .= zero(T)

      P = spzeros(T, n,n)

      q = rand(rng, T, n, 1)*T(50)
      q[end] = -one(T)
      q = vec(q)
      xtrue = rand(rng, T, n, 1)*T(50)
      strue = rand(rng, T, m, 1)*T(50)
      b = A*xtrue+strue
      b = vec(b)

      constraint = COSMO.Constraint(-A,b,COSMO.Nonnegatives)

      settings = COSMO.Settings{T}(max_iter=10000,eps_abs = 1e-5,eps_rel=1e-5)
      model = COSMO.Model{T}()
      assemble!(model,P,q,[constraint], settings = settings)
      res = COSMO.optimize!(model);

      @test res.status == :Dual_infeasible
      end
    end
  end
end
