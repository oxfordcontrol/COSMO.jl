include("../src/COSMO.jl")
using Main.COSMO, Test, LinearAlgebra, Random

rng = Random.MersenneTwister(12345)
# Conic programming example
# min 1/2*||X - S||_F
# s.t. X > 0
n = 10;
S0 = Symmetric(randn(rng, n, n))
# @testset "Lanczos Projection - Simple SDP" begin
    #for offset in [-3, 3]
        offset = 0
        S = S0 - offset*I;
        @show sum(eigen(S).values .> 0)
        q = -vec(S)
        P = SparseMatrixCSC(Diagonal(ones(n^2)));

        Aa = Diagonal(ones(n^2));
        constraints = [[
            SparseMatrixCSC(1.0*I, n^2, n^2), # I
            zeros(n^2), # 0
            COSMO.PsdCone]]

        P, q, constraints = triangulize_constraints()
        # define example problem
        settings = COSMO.Settings(verbose=false, max_iter=2500)
        model = COSMO.Model()
        COSMO.assemble!(model, P, q, constraints, settings=settings)

        res = COSMO.optimize!(model);
        E = eigen(S);
        Sp = E.vectors*max.(Diagonal(E.values),0)*E.vectors'

        @test isapprox(res.x, vec(Sp), atol=1e-3, norm=(x -> norm(x,Inf)))
    # end
# end



#=
using COSMO, Random, Test, Pkg
rng = Random.MersenneTwister(12345)
n = 1000; N = Int(n*(n + 1)/2);
P = Symmetric(sprandn(rng, n, n, 0.01))
q = randn(n)
A1 = SparseMatrixCSC(I, n, n)
b = zeros(n, n)
b = -rand(rng, m)
A = [A; -sparse(1.0I, n, n)]
b = [b; zeros(n)]

# create dual feasible problem
P = generate_pos_def_matrix(n, rng)
ytrue = rand(rng, m + n, 1)
xtrue = rand(rng, n, 1)
q = (-P * xtrue -  A' * ytrue)[:]

constraint = COSMO.Constraint(-A, b, COSMO.Nonnegatives)
=#