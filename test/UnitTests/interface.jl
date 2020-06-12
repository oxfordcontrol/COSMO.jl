using COSMO, Test, LinearAlgebra, Random, SparseArrays
rng = Random.MersenneTwister(41)

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

@testset "Interface" begin
for T in UnitTestFloats


    @testset "Interface (T = $(T))" begin
        # set!() -  to directly set the problem data without building constraints first

        P = T[4. 1;1 2]
        q = T[1; 1.]
        A = T[1. 1;1 0; 0 1]
        l = T[1.;0;0]
        u = T[1.;0.7;0.7]
        Aa = [A; -A]
        b = [u; -l]
        convex_sets = [COSMO.Nonnegatives{T}(3); COSMO.Nonnegatives{T}(3)]

        model = COSMO.Model{T}()
        COSMO.set!(model, P, q, Aa, b, convex_sets )
        res = COSMO.optimize!(model)
        @test isapprox(res.obj_val, 1.88, atol = 1e-3)

        qf = T[1; 1; 1]
        Pf = zeros(T, 1, 1)
        Af = T[1 2; 1 2]
        bf = T[1; 2]
        model = COSMO.Model{T}()
        @test_throws DimensionMismatch COSMO.set!(model, P, qf, Aa, b, convex_sets)
        @test_throws DimensionMismatch COSMO.set!(model, Pf, q, Aa, b, convex_sets)
        @test_throws DimensionMismatch COSMO.set!(model, P, q, Af, b, convex_sets)
        @test_throws DimensionMismatch COSMO.set!(model, P, q, Aa, bf, convex_sets)


        # assemble! function test different combinations of P an q

        #P not sparse, q dense
        P = rand(rng, T, 2, 2)
        q = rand(rng, T, 2)
        A = rand(rng, T, 1, 2)
        b = one(T)
        constraint = COSMO.Constraint(A, b, COSMO.Nonnegatives)
        model = COSMO.Model{T}();
        COSMO.assemble!(model, P, q, constraint)
        @test model.p.P == P
        @test model.p.q == q
        @test model.p.A == -A
        @test model.p.b == [b]

        # sparse input data
        COSMO.empty_model!(model)
        COSMO.assemble!(model, sparse(P), sparse(q), constraint)
        @test model.p.P == P
        @test model.p.q == q

        A = rand(rng, T, 5, 1)
        b = rand(rng, T, 5)
        constraint = COSMO.Constraint(A, b, COSMO.Nonnegatives)

        # P number, q vector
        COSMO.empty_model!(model)
        COSMO.assemble!(model, one(T), [one(T)], constraint)
        @test model.p.P == hcat([one(T)])
        @test model.p.q == [one(T)]

        # P vector, q number
        COSMO.empty_model!(model)
        COSMO.assemble!(model, [one(T)], one(T), constraint)
        @test model.p.P == hcat([one(T)])
        @test model.p.q == [one(T)]

        # P number, q number
        COSMO.empty_model!(model)
        COSMO.assemble!(model, one(T), one(T), constraint, settings = COSMO.Settings{T}(verbose = true))
        @test model.p.P == hcat([one(T)])
        @test model.p.q == [one(T)]

        # P vector, q matrix
        COSMO.empty_model!(model)
        COSMO.assemble!(model, [one(T)], hcat([one(T)]), constraint)
        @test model.p.P == hcat([one(T)])
        @test model.p.q == [one(T)]

        # Test case where dimension of an A in one constraint is inconsistent with P
        A = one(T)
        b = zero(T)
        P = spdiagm(0=> ones(T, 2))
        q = rand(rng, T, 2)
        constraint = COSMO.Constraint(A, b, COSMO.Nonnegatives)
        model = COSMO.Model{T}();
        @test_throws DimensionMismatch COSMO.assemble!(model, P, q, [constraint])

        # Disallow models with PSD constraints of type BigFloat
        if T == BigFloat
            model = COSMO.Model{T}();
            P = rand(rng, T, 4, 4)
            q = rand(rng, T, 4)
            A = spdiagm(0 => ones(T, 4))
            b = zeros(T, 4)
            constraint = COSMO.Constraint(A, b, COSMO.PsdCone)
            @test_throws ArgumentError COSMO.assemble!(model, P, q, [constraint])
        end
    end
end
end
