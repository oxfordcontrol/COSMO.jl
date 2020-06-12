# Primal infeasible test problems
# x >= 0, elements in A >=0, elements in b <0 and s in R+


# module PrimalInfeasibleTest
# using COSMO, Test, Random, SparseArrays, LinearAlgebra
# include("./../COSMOTestUtils.jl")

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end


@testset "Primal infeasible QP problems - Testset 1 (T = $(T))" for T in UnitTestFloats
    # for iii = 1:nn

    if precision(T) >= precision(Float64)
        nn = 1
        rng = Random.MersenneTwister(74747)

        # choose size of problem
        n = rand(rng,  5:50)
        m = 2 * n
        A = sprand(rng, T, m, n, 0.8)
        b = -rand(rng, T, m)
        A = [A; -spdiagm(0 => ones(T, n))]
        b = [b; zeros(T, n)]

        # create dual feasible problem
        P = generate_pos_def_matrix(rng, n; MT = T)
        ytrue = rand(rng, T, m + n, 1)
        xtrue = rand(rng, T, n, 1)
        q = (-P * xtrue -  A' * ytrue)[:]

        constraint = COSMO.Constraint(-A, b, COSMO.Nonnegatives)

        settings = COSMO.Settings{T}(max_iter=10000, eps_abs = 1e-5, eps_rel=1e-5)
        model = COSMO.Model{T}()
        assemble!(model, P, q, constraint, settings = settings)
        res1 = COSMO.optimize!(model);

       @test res1.status == :Primal_infeasible
   end
 end
