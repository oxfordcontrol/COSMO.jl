# Test script to test solver for a Lasso problem transformed into a second-order cone problem
workspace()
include("../Solver.jl")

using Base.Test, Convex
using OSSDP, OSSDPTypes, Mosek

# Lasso problem
# min 0.5 * ||Fz-b||_2^2 + μ||z||_1

# Generate problem data
n = 10
m = 100
F = rand(m,n)

vtrue = sprand(n,1,0.1 )
noise = 0.1*rand(m,1)
b = F*vtrue + noise
μMax = norm(F'*b,Inf)
μ = 0.1*μMax

# set up optimization problem as second order cone program of the following form
# min   0.5*w + μ*1'*t
# s.t.  -t <= z <= t
#       | 1 - w  |
#       | 2(Fz-g)|_2  <= 1 + w


#####
# Solve problem with OSSDP
#####

# create augmented matrices
Aa = [eye(n) eye(n) zeros(n,1) -eye(n) zeros(n,n) zeros(n,1) zeros(n,m+1);
    -eye(n) eye(n)  zeros(n,1) zeros(n,n) eye(n)  zeros(n,1) zeros(n,m+1);
    zeros(1,n) zeros(1,n) -1   zeros(1,n) zeros(1,n) 1 zeros(1,m+1);
    zeros(1,n) zeros(1,n) 1   zeros(1,n) zeros(1,n) 0 1 zeros(1,m);
    zeros(m,n) -2*F zeros(m,1) zeros(m,n) zeros(m,n) zeros(m,1) zeros(m,1) eye(m,m)]

ba = [zeros(n,1);zeros(n,1);1;1;-2*b]
ca = [μ*ones(n,1);zeros(n,1);0.5;zeros(2*n+1+1+m)]
Pa = zeros(length(ca),length(ca))

# define cone (t, z and w free variables, 2*n slack variables for inequality constraints, m+2 soc variables)
K = cone(2*n+1,2*n,[2+m],[])
# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0)
res = solveSDP(Pa,ca,Aa,ba,K,settings)
t = res.x[1:n]
z = res.x[n+1:n+n]
w = res.x[2n+1]


#####
# Solve problem with Convex + MOSEK
#####

t2 = Variable(n)
z2 = Variable(n)
w2 = Variable(1)
problem = minimize( 0.5*w2 + μ*sum(t2))
problem.constraints += z2 <= t2
problem.constraints += z2 >= -t2
problem.constraints += norm([1-w2;2*(F*z2-b)],2) <= 1 + w2

solve!(problem,MosekSolver())

@testset "Run all tests" begin
  @testset "Lasso (SOCP) Problem OSSDP " begin
    @test norm([1-w;2*(F*z-b)],Inf) <= 1 + w
    @test minimum(z + t) >= -1e-5
    @test maximum(z - t) <= 1e-5
  end

  @testset "Lasso (SOCP) Problem Mosek " begin
    @test norm([1-w2.value;2*(F*z2.value-b)],Inf) <= 1 + w2.value
    @test minimum(z2.value + t2.value) >= -1e-5
    @test maximum(z2.value - t2.value) <= 1e-5
  end

  @testset "Compare solvers " begin
    @test norm(z-z2.value,Inf) <= 1e-3
    @test abs(res.cost - problem.optval) <= 1e-3
  end

end
nothing

