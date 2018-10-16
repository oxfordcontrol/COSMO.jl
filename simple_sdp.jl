include("src/COSMO.jl")
using Main.COSMO
using Test, LinearAlgebra

# Conic programming example
# min 1/2*||X - S||_F
# s.t. X > 0
n = 40;
S = randn(n, n) - 4*I;
S = (S + S')/2;
Aa = Diagonal(ones(n^2));
ba = zeros(n^2,1);
c = -vec(S)
P = Diagonal(ones(n^2));

constraint1 = COSMO.Constraint(Aa,ba,COSMO.PsdCone)

# define example problem
settings = COSMO.Settings(verbose=false)

model = COSMO.Model()
assemble!(model,P,c,(constraint1))

res = COSMO.optimize!(model,settings);

E = eigen(S);
Sp = E.vectors*max.(Diagonal(E.values),0)*E.vectors'

@testset "Lanczos SDP" begin
  @test isapprox(res.x, vec(Sp), atol=1e-2, norm=(x -> norm(x,Inf)))
end

nothing