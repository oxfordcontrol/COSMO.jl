# Test of the OSSDP solver with various problems from the SDPLibrary
# SDP problem are solved in dual form (problem in primal form is defined using LMIs)
# problem: max  tr( F0*Y)
#          s.t. tr(Fi Y) = ci, for i=1,...,m
#               Y ⪴ 0

workspace()
include("../Types.jl")
include("../Scaling.jl")
include("../Projections.jl")
include("../Solver.jl")

using JLD
using Base.Test
using OSSDP, OSSDPTypes
using Convex, Mosek

data = JLD.load("./sdplib/truss1.jld")

F = data["F"]
c = data["c"]
m = data["m"]
n = data["n"]
optVal = data["optVal"]


# Solve problem with ConvexJL + Mosek
Y = Variable(n,n)
problem = maximize(trace(F[1]*Y))
for iii = 1:m
  problem.constraints += trace(F[iii+1]*Y) == c[iii]
end
constraint2 = isposdef(Y)
problem.constraints += constraint2
solve!(problem,SCSSolver(verbose=true))
Ys = Y.value



# Rewrite problem to OSSDP compatible format:
# min   1/2 x'Px + q'x
# s.t.  Ax = b
#       X ⪴ 0

P = zeros(n^2,n^2)
q = - vec(full(F[1])) # F0 is stored in F[1], minus since we want to maximize
b = c

# create A matrix: A = [vec(F1)'; vec(F2)'; ...; vec(Fm)']
A = zeros(m,n^2)
for iii = 1:m
  A[iii,:] = vec(F[iii+1])
end

settings = sdpSettings(rho=10,sigma=1.0,alpha=1.6,max_iter=10000,verbose=true)
res,dbg,sm = solveSDP(P,q,A,b,settings)
Y = reshape(res.x,n,n)

# Check if problem was solved correctly
# since we are computing -cost (min instead of max)
@testset "Compare results" begin
  @test norm(-res.cost - problem.optval) <= 1e-3
  @test maximum(abs.(Y-Ys)) < 1e-3
end

