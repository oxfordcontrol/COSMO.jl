# SOCP Lasso Testproblem

using QOCS, Test, LinearAlgebra, SparseArrays



# generate problem data
rng = Random.MersenneTwister(12345)
n = 8
m = 50*n
F = rand(rng,m,n)

vtrue = sprand(rng,n,1,0.1 )
noise = 0.1*rand(rng,m,1)
b = F*vtrue + noise
μMax = norm(F'*b,Inf)
μ = 0.1*μMax


# define lasso problem as SOCP

# create augmented matrices
Aa = sparse([1 zeros(1,2*n+1) 1 zeros(1,m);
      -1 zeros(1,2*n) 1 zeros(1,m+1);
      zeros(m,1) -2*F zeros(m,n+2) Matrix(1.0I,m,m);
      zeros(n,1) eye(n) -eye(n) zeros(n,m+2);
      zeros(n,1) -eye(n) -eye(n) zeros(n,m+2);
     zeros(1,2*n+1) -1 zeros(1,m+1);
     zeros(1,2*n+2) -1 zeros(1,m);
     zeros(m,2n+3) -eye(m)])

ba = vec([1;1;-2*b;zeros(2*n+m+2)])
q = vec([1;zeros(n);μ*ones(n,1);zeros(m+2,1)])
P = spzeros(length(q),length(q))

# define cone membership
K = QOCS.Cone(2+m,2*n,[m+2],[])

settings = QOCS.Settings()

# Solve with OSSDP
res,nothing = QOCS.solve(P,q,Aa,ba,K,settings);


@testset "SOCP - Lasso" begin
  @test res.status == :Solved
  @test isapprox(res.cost,0.4422849814458825,atol=1e-2)
end