# SOCP Lasso Testproblem

using COSMO, Test, LinearAlgebra, SparseArrays, Random

# generate problem data
rng = Random.MersenneTwister(12345)
n = 8
m = 50 * n
F = rand(rng, m, n)

vtrue = sprand(rng, n, 1, 0.1)
noise = 0.1 * rand(rng, m, 1)
b = F * vtrue + noise
μMax = norm(F' * b, Inf)
μ = 0.1 * μMax


# define lasso problem as SOCP

A1 = -sparse([1 zeros(1, 2 * n + 1) 1 zeros(1, m);
  -1 zeros(1, 2 * n) 1 zeros(1, m + 1);
  zeros(m, 1) -2 * F zeros(m, n + 2) Matrix(1.0I, m, m)])

A2 = -sparse([zeros(n, 1) Matrix(1.0I, n, n) -Matrix(1.0I, n, n) zeros(n, m + 2);
  zeros(n, 1) -Matrix(1.0I, n, n) -Matrix(1.0I, n, n) zeros(n, m + 2)])
A3 = -sparse([zeros(1, 2 * n + 1) -1 zeros(1, m + 1);
 zeros(1, 2 * n + 2) -1 zeros(1, m);
 zeros(m, 2 * n + 3) -Matrix(1.0I, m, m)])
b1 = [1; 1; -2 * b]
b2 = zeros(2 * n)
b3 = zeros(m + 2)

q = vec([1; zeros(n); μ * ones(n, 1); zeros(m + 2, 1)])
P = spzeros(length(q), length(q))

cs1 = COSMO.Constraint(A1 ,b1, COSMO.ZeroSet)
cs2 = COSMO.Constraint(A2 ,b2, COSMO.Nonnegatives)
cs3 = COSMO.Constraint(A3 ,b3, COSMO.SecondOrderCone)


# Solve with OSSDP
model = COSMO.Model()
assemble!(model, P, q, [cs1; cs2; cs3])
res = COSMO.optimize!(model);

@testset "SOCP - Lasso" begin
@test res.status == :Solved
@test isapprox(res.obj_val,0.4422849814458825,atol=1e-2)
end
nothing
