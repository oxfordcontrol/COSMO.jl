# Test script to test solver for combined sdp and socp problem
workspace()
include("../Solver.jl")

using Base.Test, Convex
using OSSDP, OSSDPTypes, Mosek

# Consider the following problem that includes a second order cone and a sdp constraint
# min c'x
# s.t. Ax = b
#       mat(x[5:13]) in S+
#       ||x[2:4]||_2 <= x[1]

# Generate problem data

A = [0. 0 0 0 1 1 1 1 1 1 1 1 1; #sum(x(5:13) = 18)
    1 0 0 -3 0 0 0 0 0 0 0 0 0; #( t = 3*x[4] )
     0  1 -1 0 0 0 0 0 0 0 0 0 0; #(x[2] = x[3])
     0 0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1;
     ]
b = [18.; 0; 0; 3; 2; 3]

c = [5.;1;1;1;3*ones(9,1)]


#####
# Solve problem with OSSDP
#####

# create augmented matrices

P = zeros(length(c),length(c))

# define cone
K = cone(0,0,[4],[9])
# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0)
res = solveSDP(P,c,A,b,K,settings)
t = res.x[1]
x1 = res.x[2:4]
x2 = res.x[5:13]


# #####
# # Solve problem with Convex + MOSEK
# #####

# t2 = Variable(n)
# z2 = Variable(n)
# w2 = Variable(1)
# problem = minimize( 0.5*w2 + μ*sum(t2))
# problem.constraints += z2 <= t2
# problem.constraints += z2 >= -t2
# problem.constraints += norm([1-w2;2*(F*z2-b)],2) <= 1 + w2

# solve!(problem,MosekSolver())

# @testset "Run all tests" begin
  @testset "Lasso (SOCP) Problem OSSDP " begin
    @test norm(x1,2) - t<= 1e-5
    @test isposdef(reshape(x2,3,3))
    @test norm(A*res.x - b,Inf) <= 1e-4
  end

#   @testset "Lasso (SOCP) Problem Mosek " begin
#     @test norm([1-w2.value;2*(F*z2.value-b)],Inf) <= 1 + w2.value
#     @test minimum(z2.value + t2.value) >= -1e-5
#     @test maximum(z2.value - t2.value) <= 1e-5
#   end

#   @testset "Compare solvers " begin
#     @test norm(z-z2.value,Inf) <= 1e-3
#     @test abs(res.cost - problem.optval) <= 1e-3
#   end

# end
# nothing

