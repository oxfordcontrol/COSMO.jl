include("src/QOCS.jl")
using Main.QOCS, Test, LinearAlgebra

# Conic programming example
# min 1/2*||X - S||_F
# s.t. X > 0
n = 4;
S = randn(n, n);
S = (S + S')/2;
Aa = Diagonal(ones(n^2));
ba = zeros(n^2,1);
c = -vec(S)
P = Diagonal(ones(n^2));

constraint1 = QOCS.Constraint(Aa,ba,QOCS.PositiveSemidefiniteCone())

# define example problem
settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,check_termination=1,eps_abs = 1e-6, eps_rel = 1e-6)

model = QOCS.Model()
assemble!(model,P,c,[constraint1])

res = QOCS.optimize!(model,settings);

E = eigen(S);
Sp = E.vectors*max.(Diagonal(E.values),0)*E.vectors'

@testset "Lanczos SDP" begin
  @test isapprox(res.x, vec(Sp), atol=1e-3, norm=(x -> norm(x,Inf)))
end

nothing
