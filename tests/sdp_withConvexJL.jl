# Test script to have an easy example solved by SCS that can be used to compare OSSDP solutions against
workspace()
include("../ossdp.jl")

using Convex, SCS
using Base.Test
using OSSDP


# # Problem DATA
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

# ############## Solution with ConvexJL = SCS
# X = Variable(3,3)

# problem = minimize(trace(C*X))
# constraint1 = trace(A1*X) ==b1
# constraint2 = trace(A2*X) ==b2
# constraint3 = lambdamin(X) >=0

# problem.constraints += [constraint1,constraint2,constraint3]

# solve!(problem,SCSSolver(verbose=false))

# # print result
# problem.status
# @show(X)



############## Solution with OSSDP
P = zeros(3,3)
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]
# define example problem
settings = sdpSettings(rho=1.0,verbose=true)

# solve SDP problem
res = solveSDP(P,q,A,b,settings)
X = reshape(res.x,3,3)

@test norm(X-[1.05593 0.369185 0.868298; 0.369185 0.129079 0.303585; 0.868299 0.303585 0.714011],Inf) < 1e-3