# Primal infeasible test problems
# Problems with m1 feasible equality constraints, one infeasible m2-sized second-order-cone constraint and one feasible psd constraint
# constrain s[m1] = t of [t;x] to be -1, which makes the second-order-cone constraint infeasible

nn = 100
rng = MersenneTwister(23123)

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
    s3 = vec(Helper.generatePosDefMatrix(r,rng))
    s = [s1;s2;s3]
    b = A*xtrue+s

    # enforce the t in ||x|| < t of second order cone constraint to be t < 0 --> infeasible
    A[m1+1,:] = zeros(n)
    b[m1+1] = -1

    # create dual feasible problem Px+q+A'y = 0, and y in K*
    P = Helper.generatePosDefMatrix(n,rng)
    ytrue_1 = randn(rng,m1,1)*50
    ytrue_2 = randn(rng,m2-1,1)*50
    ytrue_2 = [norm(ytrue_2,2)+1;ytrue_2]
    ytrue_3 = vec(Helper.generatePosDefMatrix(r,rng))
    ytrue = [ytrue_1;ytrue_2;ytrue_3]
    q = (-P*xtrue - A'*ytrue)[:]

    Kf = m1
    Kl = 0
    Kq = [m2]
    Ks = [r^2]

     K = OSSDPTypes.Cone(Kf,Kl,Kq,Ks)
     settings = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=3000,verbose=false,checkTermination=10,scaling = 10,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=true)

     res,nothing = OSSDP.solve(P,q,A,b,K,settings);
     @test res.status == :primal_infeasible
  end
end