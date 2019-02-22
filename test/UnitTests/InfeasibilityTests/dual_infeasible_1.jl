# Dual infeasible test problems
# Choose A, P and q in a way such that the problem becomes unbounded below
# here the last element of x appears in the cost function with a negative sign and is unconstrained

rng = Random.MersenneTwister(555)
nn = 1


@testset "Dual infeasible QP problems - Testset 1" begin
  for iii = 1:nn

  # choose size of problem
  n = rand(rng,5:50)
  m = 2*n
  A = sprandn(m,n,0.7)*50
  A[:,end] .= 0

  P = spzeros(n,n)

  q = randn(n,1)*50
  q[end] = -1
  q = vec(q)
  xtrue = randn(n,1)*50
  strue = rand(m,1)*50
  b = A*xtrue+strue
  b = vec(b)

  constraint = COSMO.Constraint(-A,b,COSMO.Nonnegatives)

  settings = COSMO.Settings(max_iter=10000,eps_abs = 1e-5,eps_rel=1e-5)
  model = COSMO.Model()
  assemble!(model,P,q,[constraint], settings = settings)
  res = COSMO.optimize!(model);

  @test res.status == :Dual_infeasible
  end
end