# Primal infeasible test problems
# x >= 0, elements in A >=0, elements in b <0 and s in R+

nn = 1
rng = Random.MersenneTwister(74747)

@testset "Primal infeasible QP problems - Testset 1" begin
  for iii = 1:nn

    # choose size of problem
    n = rand(rng, 5:50)
    m = 2 * n
    A = sprand(rng, m, n, 0.8)
    b = -rand(rng, m)
    A = [A; -sparse(1.0I, n, n)]
    b = [b; zeros(n)]

    # create dual feasible problem
    P = generate_pos_def_matrix(rng, n)
    ytrue = rand(rng, m + n, 1)
    xtrue = rand(rng, n, 1)
    q = (-P * xtrue -  A' * ytrue)[:]

    constraint = COSMO.Constraint(-A, b, COSMO.Nonnegatives)
    ra = 0.

    settings = COSMO.Settings(max_iter=10000, eps_abs = 1e-5, eps_rel=1e-5)
    model = COSMO.Model()
    assemble!(model, P, q, constraint, settings = settings)
    res1 = COSMO.optimize!(model);

     @test res1.status == :Primal_infeasible
  end
end
nothing



