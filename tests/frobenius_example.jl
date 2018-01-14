# Test script for a frobenius norm minimization problem
workspace()
include("../Projections.jl")
include("../Solver.jl")
include("../Helper.jl")

using Base.Test
using OSSDP
using Helper
# using PyPlot

rng = MersenneTwister(123)

# Problem data
M = rand(rng,10,10)

# Construct matrices according to solver format
n = size(M,1)
P = 2*eye(n^2)
q = -2*vec(M)
A = zeros(2,n^2)
b = [0.;0.0]


# define example problem
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.2,max_iter=500,verbose=true)
# solve SDP problem
res,dbg = solveSDP(P,q,A,b,settings)
X = reshape(res.x,n,n)

println("Is solution matrix X positive semidefinite? --> $(isNumericallyPosSemDef(X,-1e-10))")

println("||X-M||_F^2= $(vecnorm(X-M,2)^2)")