using COSMO, Test, LinearAlgebra, Random

rng = Random.MersenneTwister(12345)
# Conic programming example
# min 1/2*||X - S||_F
# s.t. X > 0
n = 10;
S0 = Symmetric(randn(rng, n, n))
@testset "Lanczos Projection - Simple SDP" begin
  for offset in [-.5, .5]
    S = S0 - offset*I;
    # @show sum(eigen(S).values .> 0)
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

    @test isapprox(res.x, vec(Sp), atol=1e-3, norm=(x -> norm(x,Inf)))
  end
end

nothing