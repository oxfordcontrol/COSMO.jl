# Dual infeasible test problems
# Problems with m1 equality constraints, m3-dimensional second order constraint and m4-dimensional psd constraint
# x1 is unbounded below --> hence problem is dual infeasible



# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

for T in UnitTestFloats
  # COSMO currently doesn't support PSD constraints with BigFloat precision
  local rng
  if precision(T) >= precision(Float64) && T != BigFloat

    nn = 1
    local rng = Random.MersenneTwister(1313)

    sum_detected = 0

    @testset "Dual infeasible QP problems - Testset 2" begin
      for iii = 1:nn

       # choose size of problem
        n = rand(rng,10:50)
        m1 = rand(rng,2:10) #equality constraints
        m2 = 1 #inequality constraints
        m3 = rand(rng,3:10) # second order cone constraint
        r = rand(rng, 4:10)
        m4 = r^2 # one sdp cone constraint of size r


        m = m1+m2+m3+m4


        # make problem primal feasible
        A = sprand(rng, T, m,n,0.8)*50
        xtrue = rand(rng,T, n,1)*50
        s1 = zeros(T, m1,1)
        s2 = rand(rng, T)
        s3 = rand(rng, T, m3-1,1)
        s3 = [norm(s3,2)+1;s3]
        s4 = vec(generate_pos_def_matrix(rng, r, MT = T))
        s = [s1;s2;s3;s4]

        # make problem unbounded in x1
        A[:,1] .= zero(T)
        A[m1+1,:] = [-one(T);zeros(n-1)]
        b = A*xtrue+s
        b[m1+1] = zero(T)
        b = vec(b)

        P = spzeros(T, n,n)
        q = vec([-one(T);rand(rng, T, n-1)])

        A1 = -A[1:m1,:]
        A2 = -A[m1+1:m1+m2,:]
        A3 = -A[m1+m2+1:m1+m2+m3,:]
        A4 = -A[m1+m2+m3+1:end,:]
        b1 = b[1:m1]
        b2 = b[m1+1:m1+m2]
        b3 = b[m1+m2+1:m1+m2+m3]
        b4 = b[m1+m2+m3+1:end]

        cs1 = COSMO.Constraint(A1,b1,COSMO.ZeroSet)
        cs2 = COSMO.Constraint(A2,b2,COSMO.Nonnegatives)
        cs3 = COSMO.Constraint(A3,b3,COSMO.SecondOrderCone)
        cs4 = COSMO.Constraint(A4,b4,COSMO.PsdCone)

        settings = COSMO.Settings{T}(max_iter=10000,eps_abs = 1e-5,eps_rel=1e-5)
        model = COSMO.Model{T}()
        assemble!(model,P,q,[cs1;cs2;cs3;cs4], settings = settings)
        res = COSMO.optimize!(model);


        @test res.status == :Dual_infeasible
      end
    end
  end
end
nothing
