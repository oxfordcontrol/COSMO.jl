using COSMO, Test, LinearAlgebra, Statistics, Random

# This test is precision agnostic and will be run with TestFloat precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end
mutable struct TestProblem
    P
    q
    constraints
end


@testset "Simple Tests" begin

for TestFloat in UnitTestFloats
    tol = TestFloat(1e-3)
    rng = Random.MersenneTwister(41)


    function simpleQP(Type::Type{T}) where {T <: AbstractFloat}
        P = T[4. 1;1 2]
        q = T[1; 1.]
        A = T[1. 1;1 0; 0 1]
        l = T[1.;0;0]
        u = T[1.;0.7;0.7]

        constraint1 = COSMO.Constraint(-A, u, COSMO.Nonnegatives)
        constraint2 = COSMO.Constraint(A, -l, COSMO.Nonnegatives)
        constraints = [constraint1; constraint2]
        return TestProblem(P, q, constraints)
    end



    @testset "Simple Tests (T = $(TestFloat))" begin

        @testset "Simple QP" begin
            p = simpleQP(TestFloat)
            model = COSMO.Model{TestFloat}()
            assemble!(model, p.P, p.q, p.constraints, settings = COSMO.Settings{TestFloat}(verbose=true))

            res = COSMO.optimize!(model);

            @test res.status == :Solved
            @test isapprox(norm(res.x - TestFloat[0.3; 0.7]), zero(TestFloat), atol=tol)
            @test isapprox(res.obj_val, TestFloat(1.8800000298331538), atol=tol)

            @test typeof(model) == COSMO.Model{TestFloat}
            @test typeof(model.p) == COSMO.ProblemData{TestFloat}
            @test typeof(model.settings) == COSMO.Settings{TestFloat}
            @test typeof(res) == COSMO.Result{TestFloat}
            @test typeof(res.x) <: AbstractVector{TestFloat}
            @test typeof(res.y) <: AbstractVector{TestFloat}
            @test typeof(res.s) <: AbstractVector{TestFloat}
        end

        @testset "update_max_iter" begin
            p = simpleQP(TestFloat)
            settings = COSMO.Settings{TestFloat}(max_iter=20)
            model = COSMO.Model{TestFloat}()
            assemble!(model, p.P, p.q, p.constraints, settings = settings)

            res = COSMO.optimize!(model);
            @test res.status == :Max_iter_reached
        end



        @testset "update_check_termination" begin
             p = simpleQP(TestFloat)
            settings = COSMO.Settings{TestFloat}(check_termination=100000)
            model = COSMO.Model{TestFloat}()
            assemble!(model, p.P, p.q, p.constraints, settings = settings)

            res = COSMO.optimize!(model);

            @test res.status == :Max_iter_reached

        end

        @testset "timelimit" begin
            p = simpleQP(TestFloat)
            settings = COSMO.Settings{TestFloat}(time_limit = 1, check_termination = 100000000, max_iter = 10000000)
            model = COSMO.Model{TestFloat}()
            assemble!(model, p.P, p.q, p.constraints, settings = settings)

            res = COSMO.optimize!(model);
            @test res.status == :Time_limit_reached
        end

         @testset "warm_starting" begin
            p = simpleQP(TestFloat)
            settings = COSMO.Settings{TestFloat}(check_termination = 1)
            model = COSMO.Model{TestFloat}()
            assemble!(model, p.P, p.q, p.constraints, settings = settings)

            res1 = COSMO.optimize!(model);
            n = 2
            m = 6
            x0 = res1.x
            y0 = res1.y
            s0 = res1.s
            model = COSMO.Model{TestFloat}()
            assemble!(model, p.P, p.q, p.constraints, settings = settings)
            COSMO.warm_start_primal!(model, x0)
            COSMO.warm_start_dual!(model, y0)

            # check that warm starting was applied correctly
            @test model.vars.x == x0
            @test model.vars.μ == -y0
            @test norm(model.vars.s.data - s0, 2) < 1e-3

            # actually, warm start with a bit of noise
            COSMO.warm_start_primal!(model, x0 .+ TestFloat(0.01) * rand(rng, TestFloat, n))
            COSMO.warm_start_dual!(model, y0 .+ TestFloat(0.01) * rand(rng, TestFloat, m))

            res2 = COSMO.optimize!(model);

            @test res1.status == :Solved && res2.status == :Solved && res2.iter < res1.iter

            @test_throws DimensionMismatch COSMO.warm_start_primal!(model, rand(TestFloat,4))
            @test_throws DimensionMismatch COSMO.warm_start_dual!(model, rand(TestFloat, 2))
        end

    end
end # UnitTestFloats
end
nothing
