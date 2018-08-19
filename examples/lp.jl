# Test script to test solver for an lp

using Test
using QOCS, LinearAlgebra

# Linear program example
# min c'x
# s.t. Ax <= b
#      x >= 1,  x2 =>5, x1+x3 => 4
c = [1; 2; 3; 4]
A = Matrix(1.0I,4,4)
b = [10; 10; 10; 10]

# create augmented matrices
Aa = -[A;-Matrix(1.0I,4,4);0 -1 0 0;-1 0 -1 0]
ba = vec([b; -ones(4,1);-5;-4])
P = zeros(size(A,2),size(A,2))

constraint1 = QOCS.Constraint(Aa,ba,QOCS.Nonnegatives())

# define example problem
settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,check_termination=1,eps_abs = 1e-6, eps_rel = 1e-6)

model = QOCS.Model()
assemble!(model,P,c,[constraint1])

res, = QOCS.optimize!(model,settings);

@testset "Linear Problem" begin
  @test isapprox(res.x[1:4],[3;5;1;1], atol=1e-2, norm=(x -> norm(x,Inf)))
  @test isapprox(res.cost,20.0, atol=1e-2)
end
nothing