# Test script to test solver for an lp
workspace()
include("../src/Solver.jl")

using Base.Test
using OSSDP, OSSDPTypes

# Linear program example
# min c'x
# s.t. Ax <= b
#      x >= 1,  x2 =>5, x1+x3 => 4
c = [1; 2; 3; 4]
A = eye(4)
b = [10; 10; 10; 10]

# create augmented matrices
Aa = [A eye(4) zeros(4,6);eye(4) zeros(4,4) -eye(4) zeros(4,2); 0 1 0 0 zeros(1,8) -1 0; 1 0 1 0 zeros(1,8) 0 -1]
ba = [b; ones(4,1);5;4]
ca = [c;zeros(10,1)]
P = zeros(size(Aa,2),size(Aa,2))
# define cone
K = OSSDPTypes.Cone(4,10,[],[])
# define example problem
settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0)
res,  = OSSDP.solve(P,ca,Aa,ba,K,settings)

@testset "Linear Problem" begin
  @test isapprox(res.x[1:4],[3;5;1;1], atol=1e-4, norm=(x -> norm(x,Inf)))
  @test isapprox(res.cost,20.0, atol=1e-4)
end
