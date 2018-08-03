# QP Lasso testproblem

using OSSDP, Base.Test



# generate problem data
rng = MersenneTwister(1313)
n = 15
m = 50*n
A = sprandn(rng,m,n,0.5)
vtrue = 1/n*sprand(rng,n,0.5)
noise = 1/4*randn(rng,m)
b = full(A*vtrue + noise)
λ = 0.2*norm(A'*b,Inf)


# define lasso problem as QP
Aa = [-A zeros(m,n) eye(m,m);
       eye(n,n) -eye(n,n) zeros(n,m);
       -eye(n,n) -eye(n,n) zeros(n,m)]

ba = [-b;zeros(2*n)]
P = 2*diagm([zeros(2*n);ones(m)]) # times two to cancel the 1/2 in the cost function
q = [zeros(n);λ*ones(n);zeros(m)]
K = OSSDPTypes.Cone(m,2*n,[],[])

settings = OSSDPTypes.OSSDPSettings()

res,nothing = OSSDP.solve(P,q,Aa,ba,K,settings);


@testset "QP - Lasso" begin
  @test res.status == :solved
  @test isapprox(res.cost, 46.40521553063313, atol=1e-3)
end
