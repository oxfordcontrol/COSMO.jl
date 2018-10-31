# QP Lasso testproblem

using Test, LinearAlgebra, SparseArrays, Random
using COSMO

# generate problem data
rng = Random.MersenneTwister(1313)
n = 15
m = 50*n

A = sprandn(rng,m,n,0.5)
vtrue = 1/n*sprand(rng,n,0.5)
noise = 1/4*randn(rng,m)
b = Vector(A*vtrue + noise)
λ = 0.2*norm(A'*b,Inf)


# define lasso problem as QP
A1 = [A zeros(m,n) -Matrix(1.0I,m,m)]
A2 = [-Matrix(1.0I,n,n) Matrix(1.0I,n,n) zeros(n,m);
       Matrix(1.0I,n,n) Matrix(1.0I,n,n) zeros(n,m)]

b1 = -b;
b2 = zeros(2*n)


P = 2*Matrix(Diagonal([zeros(2*n);ones(m)]))# times two to cancel the 1/2 in the objVal function
q = [zeros(n);λ*ones(n);zeros(m)]


constraint1 = COSMO.Constraint(A1,b1,COSMO.ZeroSet)
constraint2 = COSMO.Constraint(A2,b2,COSMO.Nonnegatives)
constraints = [constraint1;constraint2]

settings = COSMO.Settings()
model = COSMO.Model()
assemble!(model,P,q,constraints)

res = COSMO.optimize!(model,settings)

@testset "QP - Lasso" begin
  @test res.status == :Solved
  @test isapprox(res.objVal, 46.40521553063313, atol=1e-3)
end
nothing
