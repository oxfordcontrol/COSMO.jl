
using COSMO, Test, LinearAlgebra, Statistics, Random

tol = 1e-3
rng = Random.MersenneTwister(41)
mutable struct TestProblem
    P
    q
    constraints
end

function simpleQP()
    P = [4. 1;1 2]
    q = [1; 1.]
    A = [1. 1;1 0; 0 1]
    l = [1.;0;0]
    u = [1.;0.7;0.7]

    constraint1 = COSMO.Constraint(-A, u, COSMO.Nonnegatives)
    constraint2 = COSMO.Constraint(A, -l, COSMO.Nonnegatives)
    constraints = [constraint1; constraint2]
    return TestProblem(P, q, constraints)
end



@testset "Simple Tests" begin

    @testset "Simple QP" begin
        p = simpleQP()
        model = COSMO.Model()
        assemble!(model, p.P, p.q, p.constraints, settings = COSMO.Settings(verbose=true))

        res = COSMO.optimize!(model);

        @test res.status == :Solved
        @test isapprox(norm(res.x - [0.3; 0.7]), 0., atol=tol)
        @test isapprox(res.obj_val, 1.8800000298331538, atol=tol)

    end

    @testset "update_max_iter" begin
        p = simpleQP()
        settings = COSMO.Settings(max_iter=20)
        model = COSMO.Model()
        assemble!(model, p.P, p.q, p.constraints, settings = settings)

        res = COSMO.optimize!(model);
        @test res.status == :Max_iter_reached
    end



    @testset "update_check_termination" begin
         p = simpleQP()
        settings = COSMO.Settings(check_termination=100000)
        model = COSMO.Model()
        assemble!(model, p.P, p.q, p.constraints, settings = settings)

        res = COSMO.optimize!(model);

        @test res.status == :Max_iter_reached

    end


    @testset "update_rho" begin

        # TODO: write problem where rho is updated, results should be <exactly> the same

    end


    @testset "timelimit" begin
        p = simpleQP()
        settings = COSMO.Settings(time_limit=1, check_termination=100000000,max_iter=10000000)
        model = COSMO.Model()
        assemble!(model, p.P, p.q, p.constraints, settings = settings)

        res = COSMO.optimize!(model);
        @test res.status == :Time_limit_reached
    end

     @testset "warm_starting" begin
        p = simpleQP()
        settings = COSMO.Settings(check_termination = 1)
        model = COSMO.Model()
        assemble!(model, p.P, p.q, p.constraints, settings = settings)

        res1 = COSMO.optimize!(model);
        n = 2
        m = 6
        x0 = res1.x
        y0 = res1.y
        s0 = res1.s

        model = COSMO.Model()
        assemble!(model, p.P, p.q, p.constraints, settings = settings)
        COSMO.warm_start_primal!(model, x0)
        COSMO.warm_start_dual!(model, y0)
        COSMO.warm_start_slack!(model, s0)
        res2 = COSMO.optimize!(model);

        @test res1.status == :Solved && res2.status == :Solved && res2.iter < res1.iter

        @test_throws DimensionMismatch COSMO.warm_start_primal!(model, rand(4))
        @test_throws DimensionMismatch COSMO.warm_start_dual!(model, rand(2))
    end

end
nothing
