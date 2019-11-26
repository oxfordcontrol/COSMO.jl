# Primal infeasible test problems
# Problems with m1 feasible equality constraints, one infeasible m2-sized second-order-cone constraint and one feasible psd constraint
# constrain s[m1] = t of [t;x] to be -1, which makes the second-order-cone constraint infeasible

nn = 1
rng = Random.MersenneTwister(1313)

@testset "Primal infeasible QP problems - Testset 3" begin
  for iii = 1:nn

    # choose size of problem
    n = rand(rng,10:50)
    m1 = rand(rng,2:10) #equality constraint
    m2 = rand(rng,3:10) # second order cone constraint
    r = rand(rng,4:10)
    m3 = r^2 # one sdp cone constraint of size r
    m = m1+m2+m3

    A = sprandn(rng,m,n,0.8)*50
    xtrue = rand(rng,n,1)*50
    s1 = zeros(m1,1)
    s2 = randn(rng,m2,1)
    s3 = vec(generate_pos_def_matrix(rng, r))
    s = [s1;s2;s3]
    b = A*xtrue+s

    # enforce the t in ||x|| < t of second order cone constraint to be t < 0 --> infeasible
    A[m1+1,:] = zeros(n)
    b[m1+1] = -1
    b = vec(b)

    # create dual feasible problem Px+q+A'y = 0, and y in K*
    P = generate_pos_def_matrix(rng, n)
    ytrue_1 = randn(rng,m1,1)*50
    ytrue_2 = randn(rng,m2-1,1)*50
    ytrue_2 = [norm(ytrue_2,2)+1;ytrue_2]
    ytrue_3 = vec(generate_pos_def_matrix(rng, r))
    ytrue = [ytrue_1;ytrue_2;ytrue_3]
    q = (-P*xtrue - A'*ytrue)[:]

    A1 = -A[1:m1,:]
    A2 = -A[m1+1:m1+m2,:]
    A3 = -A[m1+m2+1:end,:]
    b1 = b[1:m1,:]
    b2 = b[m1+1:m1+m2]
    b3 = b[m1+m2+1:end]

    cs1 = COSMO.Constraint(A1,b1,COSMO.ZeroSet)
    cs2 = COSMO.Constraint(A2,b2,COSMO.SecondOrderCone)
    cs3 = COSMO.Constraint(A3,b3,COSMO.PsdCone)

    settings = COSMO.Settings(max_iter=10000,eps_abs = 1e-5,eps_rel=1e-5)

    model = COSMO.Model()
    assemble!(model,P,q,[cs1;cs2;cs3], settings = settings)
    res = COSMO.optimize!(model);

    @test res.status == :Primal_infeasible
  end
end
nothing
