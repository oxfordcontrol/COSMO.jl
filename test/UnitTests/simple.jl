
using QOCS, Test, LinearAlgebra, Statistics, Random

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

    constraint1 = QOCS.Constraint(-A,u,QOCS.Nonnegatives())
    constraint2 = QOCS.Constraint(A,-l,QOCS.Nonnegatives())
    constraints = [constraint1;constraint2]
    return TestProblem(P,q,constraints)
end



@testset "Simple Tests" begin

    @testset "Simple QP" begin
        p = simpleQP()
        settings = QOCS.Settings()
        model = QOCS.Model()
        assemble!(model,p.P,p.q,p.constraints)

        res = QOCS.optimize!(model,settings);

        @test res.status == :Solved
        @test isapprox(norm(res.x - [0.3; 0.7]), 0., atol=tol)
        @test isapprox(res.objVal, 1.8800000298331538, atol=tol)

    end


    # @testset "Update_b" begin
    #     p = simpleQP()
    #     p.b = p.b*0.9
    #     settings = QOCS.Settings()
    #     res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

    #     @test res.status == :Solved
    #     @test isapprox(norm(res.x - [0.27; 0.63]), 0., atol=tol)
    #     @test isapprox(res.objVal, 1.6128000168085233, atol=tol)

    # end

    # @testset "Update_q" begin
    #     p = simpleQP()
    #     p.q = [-10;10]
    #     settings = QOCS.Settings()
    #     res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

    #     @test res.status == :Solved
    #     @test isapprox(norm(res.x - [0.7; 0.3]), 0., atol=tol)
    #     @test isapprox(res.objVal, -2.7199998274697608, atol=tol)

    # end

    @testset "update_max_iter" begin
        p = simpleQP()
        settings = QOCS.Settings(max_iter=20)
        model = QOCS.Model()
        assemble!(model,p.P,p.q,p.constraints)

        res = QOCS.optimize!(model,settings);


        @test res.status == :Max_iter_reached
    end



    @testset "update_check_termination" begin
         p = simpleQP()
        settings = QOCS.Settings(check_termination=100000)
        model = QOCS.Model()
        assemble!(model,p.P,p.q,p.constraints)

        res = QOCS.optimize!(model,settings);

        @test res.status == :Max_iter_reached

    end


    @testset "update_rho" begin

        # TODO: write problem where rho is updated, results should be <exactly> the same

    end


    @testset "timelimit" begin
        p = simpleQP()
        settings = QOCS.Settings(time_limit=1, check_termination=100000000,max_iter=10000000)
        model = QOCS.Model()
        assemble!(model,p.P,p.q,p.constraints)

        res = QOCS.optimize!(model,settings);
        @test res.status == :Time_limit_reached
    end

     @testset "warm_starting" begin
        p = simpleQP()
        settings = QOCS.Settings(check_termination = 1)
        model = QOCS.Model()
        assemble!(model,p.P,p.q,p.constraints)

        res1 = QOCS.optimize!(model,settings);
        n = 2
        m = 6
        x0 = res1.x + 0.01*randn(rng,n)
        y0 = res1.y + 0.01*randn(rng,m)
        warmStart!(model,x0=x0,y0=y0)
        res2 = QOCS.optimize!(model,settings);

        @test res1.status == :Solved && res2.status == :Solved && res2.iter < res1.iter
    end

end
nothing