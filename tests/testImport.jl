# Quick test to make sure that sdpLibImport.jl works as intended
# Solve problem defined by imported file using ConvexJL and SCS
workspace()
using Convex, SCS, JLD
using Base.Test


data = JLD.load("./sdplib/truss2.jld")

F = data["F"]
c = data["c"]
m = data["m"]
n = data["n"]
optVal = data["optVal"]


# Describe problem using SCS
x = Variable(m)
X = Variable(n,n)
problem = minimize(c'*x)
constraint1 = sum(F[i+1]*x[i] for i=1:m) - F[1] == X
constraint2 = isposdef(X)
problem.constraints += [constraint1,constraint2]
solve!(problem,SCSSolver(verbose=false))

@test norm(problem.optval - optVal) <= 1e-3
