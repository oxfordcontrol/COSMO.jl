
using QOCS, Base.Test
tol = 1e-3

mutable struct TestProblem
    P
    q
    A
    b
    K
end

function simpleQP()
    P = [4. 1;1 2]
    q = [1; 1.]
    A = [1. 1;1 0; 0 1]
    l = [1.;0;0]
    u = [1.;0.7;0.7]
    Aa = [A;-A]
    b = [u;-l]
    K = QOCS.Cone(0,6,[],[])
    return TestProblem(P,q,Aa,b,K)
end



@testset "Simple Tests" begin

    @testset "Simple QP" begin
        p = simpleQP()
        settings = QOCS.Settings()
        res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

        @test res.status == :Solved
        @test isapprox(norm(res.x - [0.3; 0.7]), 0., atol=tol)
        @test isapprox(res.cost, 1.8800000298331538, atol=tol)

    end


    @testset "Update_b" begin
        p = simpleQP()
        p.b = p.b*0.9
        settings = QOCS.Settings()
        res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

        @test res.status == :Solved
        @test isapprox(norm(res.x - [0.27; 0.63]), 0., atol=tol)
        @test isapprox(res.cost, 1.6128000168085233, atol=tol)

    end

    @testset "Update_q" begin
        p = simpleQP()
        p.q = [-10;10]
        settings = QOCS.Settings()
        res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

        @test res.status == :Solved
        @test isapprox(norm(res.x - [0.7; 0.3]), 0., atol=tol)
        @test isapprox(res.cost, -2.7199998274697608, atol=tol)

    end

    @testset "update_max_iter" begin
        p = simpleQP()
        settings = QOCS.Settings(max_iter=20)
        res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

        @test res.status == :Max_iter_reached
    end



    @testset "update_check_termination" begin
        p = simpleQP()
        settings = QOCS.Settings(check_termination=100000)
        res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

        @test res.status == :Max_iter_reached

    end


    @testset "update_rho" begin

        # TODO: write problem where rho is updated, results should be <exactly> the same

    end


    @testset "timelimit" begin
        p = simpleQP()
        settings = QOCS.Settings(timelimit=1, check_termination=100000000,max_iter=10000000)
        res,nothing = QOCS.solve(p.P,p.q,p.A,p.b,p.K,settings);

        @test res.status == :Time_limit_reached
    end

end