        # Exponential cone constrained problems

using Test, LinearAlgebra, SparseArrays
using COSMO

@testset "Exponential cone problems" begin

    @testset "Feasible" begin

        # solve the following exponential cone problem
        # max  x
        # s.t. y * exp(x / y) <= z
        #      y == 1, z == exp(5)

        P = spzeros(3, 3)
        q = [-1.; 0; 0]

        # y * exp( x / y) <= z
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = zeros(3)

        # y == 1 and z == exp(5) = 148.4131591025766
        A2 = sparse([0 1. 0; 0 0 1])
        b2 = [-1.; -exp(5)]
        cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)
        cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)

        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2])

        res = COSMO.optimize!(model)
        @test res.status == :Solved
        @test isapprox(res.obj_val, -5, atol=1e-3)
    end


    @testset "Primal infeasible 1" begin
        # min  x
        # s.t. y * exp( x / y) <= z
        #       y == 1
        #       z == -1

        P = spzeros(3, 3)
        q = [1.; 0; 0]

        # y * exp( x / y) <= z
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = zeros(3)

        # y == 1
        A2 = [0 -1. 0]
        b2 = [-1.]

        # z == -1
        A3 = [0 0 -1.]
        b3 = [1]

        cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)
        cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)
        cs3 = COSMO.Constraint(A3, b3, COSMO.ZeroSet)

        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2; cs3], settings = COSMO.Settings())

        res = COSMO.optimize!(model)
        @test res.status == :Primal_infeasible
    end

    @testset "Primal infeasible 2" begin
        # min  x
        # s.t. y * exp( x / y) <= z - 0.2
        #      y * exp( x / y) >= z + 0.2

        P = spzeros(3, 3)
        q = [1.; 0; 0]

        # y * exp( x / y) <= z - 0.2
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = [0.; 0; -0.2]

        # y * exp( x / y) >= z + 0.3
        A2 = sparse(-Matrix(1.0I, 3, 3))
        b2 = [0.; 0; -0.3]

        cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)
        cs2 = COSMO.Constraint(A2, b2, COSMO.ExponentialCone)


        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2], settings = COSMO.Settings())

        res = COSMO.optimize!(model)
        @test res.status == :Primal_infeasible
    end

    @testset "Dual infeasible" begin
        # max z
        # s.t. y * exp(x / y) <= z

        P = spzeros(3, 3)
        q = [0; 0; -1.0]

        # y * exp( x / y) <= z
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = zeros(3)

        cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)

        model = COSMO.Model()
        assemble!(model, P, q, cs1)

        res = COSMO.optimize!(model)
        @test res.status == :Dual_infeasible
    end

end

@testset "Dual Exponential cone problems" begin

    @testset "Feasible" begin

        # solve the following dual exponential cone problem
        # min  y
        # s.t. -x * exp(y / x) <= z
        #      x == -1, z == exp(5)

        P = spzeros(3, 3)
        q = [0; 1.; 0]

        # -x * exp( y / x) <= z
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = zeros(3)
        cs1 = COSMO.Constraint(A1, b1, COSMO.DualExponentialCone)

        # y == 1 and z == exp(5)
        A2 = sparse([1. 0 0; 0 0 1])
        b2 = [1.; -exp(5)]
        cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)

        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2])

        res = COSMO.optimize!(model)
        @test res.status == :Solved
        @test isapprox(res.obj_val, -6., atol=1e-3)
    end
end


nothing
