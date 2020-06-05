using COSMO, Test, LinearAlgebra, Statistics, Random

# This test is precision agnostic and will be run with TestFloat precision
if @isdefined UnitTestFloat
    TestFloat = UnitTestFloat
else
    TestFloat = Float64
end
TestFloat = BigFloat
tol = TestFloat(1e-3)
rng = Random.MersenneTwister(41)

mutable struct TestProblem
    P
    q
    constraints
end

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



@testset "Simple Tests" begin

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
        COSMO.warm_start_slack!(model, s0)
        res2 = COSMO.optimize!(model);

        @test res1.status == :Solved && res2.status == :Solved && res2.iter < res1.iter

        @test_throws DimensionMismatch COSMO.warm_start_primal!(model, rand(TestFloat,4))
        @test_throws DimensionMismatch COSMO.warm_start_dual!(model, rand(TestFloat, 2))
    end

end
nothing
