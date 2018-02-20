# Test script to test solver for combined sdp and socp problem
workspace()
include("../src/Solver.jl")
include("../src/Helper.jl")

using Base.Test, Convex, Helper
using OSSDP, OSSDPTypes, Mosek, SCS

# Consider the following problem that includes a second order cone and a sdp constraint
# min c'x
# s.t. Ax = b
#       mat(x[5:13]) in S+
#       ||x[2:4]||_2 <= x[1]

# Generate problem data

A = [0 1 0 0 0 0 0 0 0 0 0 0 0; #( x[2] = 2 )
     0 2 -1 0 0 0 0 0 0 0 0 0 0; #(2x[2] = x[3])
     0 1 0 0 -1 0 0 0 0 0 0 0 0; #x[2] = X[1,1]
     0 0 1 0 0 0 0 0 -1 0 0 0 0; #x[3] = X[2,2]
     0 0 0 1 0 0 0 0 0 0 0 0 -1; #x[4] = X[3,3]
     ]
b = [2;0;0;0;0]

c = [5.;1;1;0;zeros(9,1)]


#####
# Solve problem with OSSDP
#####

# create augmented matrices

P = zeros(length(c),length(c))

# define cone
K = Cone(0,0,[4],[9])
# define example problem
settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0)
res, = OSSDP.solve(P,c,A,b,K,settings)
t = res.x[1]
x1 = res.x[2:4]
x2 = reshape(res.x[5:13],3,3)


# #####
# # Solve problem with Convex + MOSEK
# #####

x = Variable(4)
X = Variable(3,3)

problem = minimize( c[1:4]'*x)
problem.constraints += x[2] == 2
problem.constraints += 2*x[2] == x[3]
problem.constraints += X[1,1] == x[2]
problem.constraints += X[2,2] == x[3]
problem.constraints += X[3,3] == x[4]

problem.constraints += norm(x[2:4],2) <= x[1]
problem.constraints +=  isposdef(X)


solve!(problem,MosekSolver())

 @testset "Run all tests" begin
  @testset "Combined SOCP SDP Problem - OSSDP" begin
    @test abs(norm(x1,2) - t)<= 1e-5
    @test isNumericallyPosSemDef(reshape(x2,3,3),1e-6)
    @test norm(A*res.x - b,Inf) <= 1e-4
  end

  @testset "Combined SOCP SDP Problem - Mosek " begin
    @test abs(norm(x.value[2:4],2) - x.value[1]) <= 1e-5
    @test isNumericallyPosSemDef(X.value,1e-6)
  end

  @testset "Compare solvers " begin
    @test norm(res.x[1:4]-x.value,Inf) <= 1e-3
    @test abs(res.cost - problem.optval) <= 1e-3
  end

end
# nothing

