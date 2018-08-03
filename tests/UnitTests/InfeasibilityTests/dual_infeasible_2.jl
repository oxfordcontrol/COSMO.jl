# Dual infeasible test problems
# Problems with m1 equality constraints, m3-dimensional second order constraint and m4-dimensional psd constraint
# x1 is unbounded below --> hence problem is dual infeasible


nn = 1
rng = MersenneTwister(1313)

sum_detected = 0
@testset "Dual infeasible QP problems - Testset 2" begin
  for iii = 1:nn

   # choose size of problem
    n = rand(rng,10:50)
    m1 = rand(rng,2:10) #equality constraints
    m2 = 1 #inequality constraints
    m3 = rand(rng,3:10) # second order cone constraint
    r = rand(rng,4:10)
    m4 = r^2 # one sdp cone constraint of size r


    m = m1+m2+m3+m4


    # make problem primal feasible
    A = sprandn(rng,m,n,0.8)*50
    xtrue = rand(rng,n,1)*50
    s1 = zeros(m1,1)
    s2 = rand(rng)
    s3 = randn(rng,m3-1,1)
    s3 = [norm(s3,2)+1;s3]
    s4 = vec(Helper.generatePosDefMatrix(r,rng))
    s = [s1;s2;s3;s4]

    # make problem unbounded in x1
    A[:,1] = 0
    A[m1+1,:] = [-1;zeros(n-1)]
    b = A*xtrue+s
    b[m1+1] = 0
    b = vec(b)

    P = spzeros(n,n)
    q = vec([-1;randn(rng,n-1)])


    Kf = m1
    Kl = m2
    Kq = [m3]
    Ks = [r^2]

   K = QOCS.Cone(Kf,Kl,Kq,Ks)
   setOFF = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,check_termination=1,scaling = 10,eps_prim_inf=1e-4,eps_dual_inf=1e-4,adaptive_rho=true)
   res,ws,δx,δμ = QOCS.solve(P,q,A,b,K,setOFF);

    @test res.status == :Dual_infeasible
  end
end

