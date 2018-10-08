# Test routine to test the infeasibility example from the OSQP Paper

# workspace()
# include("../../../src/Solver.jl")

using COSMO, Base.Test
maxIter = 3000
rng = Random.MersenneTwister(1313)

setUnNonAdaptive =   Settings(max_iter=maxIter,verbose=true,check_termination=1,checkInfeasibility=50,scaling = 0 ,scaleFunc=1,adaptive_rho=false,eps_abs=1e-7,eps_rel=1e-7)
setMeanNonAdaptive = Settings(max_iter=maxIter,verbose=false,check_termination=1,checkInfeasibility=50,scaling = 10,scaleFunc=2,adaptive_rho=false,eps_abs=1e-7,eps_rel=1e-7)

setUnAdaptive =   Settings(max_iter=maxIter,verbose=false,check_termination=1,checkInfeasibility=50,scaling = 0 ,scaleFunc=1,adaptive_rho=true,eps_abs=1e-7,eps_rel=1e-7)
setMeanAdaptive = Settings(max_iter=maxIter,verbose=false,check_termination=1,checkInfeasibility=50,scaling = 10,scaleFunc=2,adaptive_rho=true,eps_abs=1e-7,eps_rel=1e-7)

P = sparse([1. 0;0 0])
q = [1.;-1]

# # optimality
# a = 1
# u1 = 5
# u3 = 3
# A = sparse([1 a;-1 -a;1 0;-1 0; 0 1; 0 -1])
# b = [u1;0;3;-1;u3;-1]
# P = sparse([1 0;0 0])
# q = [1;-1]
# Kf = 0
# Kl = 6


# primal infeasibility
# case = "PRIMAL"
# A = sparse([1. 1;-1 -1;1 0;-1 0;0 1;0 -1])
# b = [0.;0;3;-1;3;-1]
# Kf = 0
# Kl = 6

# # dual infeasibility
# case = "DUAL"
# A = sparse([1 0;-1 0;1 0;-1 0; 0 -1])
# b = [2;0;3;-1;-1]
# Kf = 0
# Kl = 5

# primal and dual infeasiblity
case = "BOTH"
A = sparse([1 1;-1 -1;1 0;-1 0;;0 -1])
b = [0;0;3;-1;-1]
Kf = 0
Kl = 5


## primal infeasible test problem from MOSEK doc, i.e. min x1 s.t. x1 <= 5
# P = zeros(2,2)
# q = [1;0]
# A = [0 1;1 0]
# b = [0;5]
# Kf = 1
# Kl = 1
# Ks = []
# Kq = []


Kq = []
Ks = []
# define cone membership
K = COSMO.Cone(Kf,Kl,Kq,Ks)


# Solve with COSMO
res1,ws,δx1,δμ1 = COSMO.solve(P,q,A,b,K,setMeanAdaptive);
print("\n. δx1: $(δx1), δμ1: $(δμ1)")
# res2,nothing,δx2,δμ2 = COSMO.solve(P,q,A,b,K,setMeanNonAdaptive);
# print("\n. δx2: $(δx2), δμ2: $(δμ2)")
# res3,nothing,δx3,δμ3 = COSMO.solve(P,q,A,b,K,setUnAdaptive);
# print("\n. δx3: $(δx3), δμ3: $(δμ3)")
# res4,nothing,δx4,δμ4 = COSMO.solve(P,q,A,b,K,setMeanAdaptive);
# print("\n. δx4: $(δx4), δμ4: $(δμ4)\n")

# eps_pinf = 1e-4
# eps_dinf = 1e-4

# @testset "Test infeasibility conditions" begin
#     if case == "PRIMAL"
#       δy = -δμ1
#       @test norm(A'*δy,Inf)/norm(δy,Inf) <= eps_pinf
#       support_δy = b'*δy/norm(δy,Inf)
#       @test support_δy <= -eps_pinf
#     elseif case == "DUAL"
#       @test norm(P*δx1,Inf)/norm(δx1,Inf) <= eps_dinf
#       @test maximum(A*δx1)/norm(δx1,Inf) <= eps_dinf
#       @test q'*δx1/norm(δx1,Inf) <= -eps_dinf
#     elseif case == "BOTH"
#       @test norm(P*δx1,Inf)/norm(δx1,Inf) <= eps_dinf
#       @test maximum(A*δx1)/norm(δx1,Inf) <= eps_dinf
#       @test q'*δx1/norm(δx1,Inf) <= -eps_dinf
#       @test norm(A'*δy,Inf)/norm(δy,Inf) <= eps_pinf
#       support_δy = b'*δy/norm(δy,Inf)
#       @test support_δy <= -eps_pinf
#     end
# end


