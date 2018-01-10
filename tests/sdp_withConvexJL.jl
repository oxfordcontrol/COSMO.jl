# Test script to have an easy example solved by SCS that can be used to compare OSSDP solutions against
workspace()
include("../ossdp.jl")

using Convex, SCS
using Base.Test
using OSSDP
# using PyPlot


# # Problem DATA
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

############## Solution with ConvexJL = SCS
X = Variable(3,3)

problem = minimize(trace(C*X))
constraint1 = trace(A1*X) ==b1
constraint2 = trace(A2*X) ==b2
constraint3 = lambdamin(X) >=0

problem.constraints += [constraint1,constraint2,constraint3]

solve!(problem,SCSSolver(verbose=true))


# print result
# problem.status
# @show(problem.optval)
# @show(X.value)

# A = [vec(A1)';vec(A2)']
# b = [b1;b2]
# @test A*vec(x) == b


############## Solution with OSSDP
P = zeros(3,3)
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]
# define example problem
settings = sdpSettings(rho=1.0,sigma=100.0,alpha=1.6,max_iter=2500,verbose=true)

# solve SDP problem
res,dbg = solveSDP(P,q,A,b,settings)
X = reshape(res.x,3,3)

#plot debugging variables

# PyPlot.figure(1,facecolor="white")
# PyPlot.plot(1:1:100,dbg.ν[1:100],label="x", linewidth=2.0)
# PyPlot.plot(1:1:100,dbg.μ[1:100],label="x", linewidth=2.0)
# PyPlot.hold(true)
# PyPlot.grid(true)
# PyPlot.legend()
# PyPlot.xlabel(L"i")
# PyPlot.ylabel(L"var")
# PyPlot.title("Primal and dual variables over time")
# PyPlot.tight_layout()

# compute smallest eigenvalue of result X
F = eigfact(X)
XλMin = minimum(F[:values])
F = eigfact(reshape(res.s,3,3))
SλMin = minimum(F[:values])
epsEV = -1e-15
@testset "Compare results" begin
  @test norm(X-[1.05593 0.369185 0.868298; 0.369185 0.129079 0.303585; 0.868299 0.303585 0.714011],Inf) < 1e-3
  @test norm(res.cost-13.902255839911007) < 1e-3
  @test norm(res.x-res.s,Inf)<1e-3
  @test XλMin > epsEV
  @test SλMin > epsEV

end


