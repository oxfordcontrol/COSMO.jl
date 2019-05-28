# Power cone constrained problems
using Test, LinearAlgebra, SparseArrays
using COSMO

@testset "Power cone problems" begin

    @testset "Feasible" begin

        # solve the following power cone problem
        # max  x1^0.6 y^0.4 + x2^0.1
        # s.t. x1, y, x2 >= 0
        #      x1 + 2y  + 3x2 == 3
        # which is equivalent to
        # max z1 + z2
        # s.t. (x1, y, z1) in K_pow(0.6)
        #      (x2, 1, z2) in K_pow(0.1)
        #      x1 + 2y + 3x2 == 3

        # x = (x1, y, z1, x2, y2, z2)
        n = 6
        P = spzeros(n, n)
        q = zeros(6)
        q[3] = q[6] = -1.

        # (x1, y, z1) in K_pow(0.6)
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = zeros(3)
        cs1 = COSMO.Constraint(A1, b1, COSMO.PowerCone(0.6), 6, 1:3)

        # (x2, 2, z2) in K_pow(0.1)
        A2 = sparse(Matrix(1.0I, 3, 3))
        b2 = zeros(3)
        cs2 = COSMO.Constraint(A2, b2, COSMO.PowerCone(0.1),6, 4:6)

        # x1 + 2y + 3x2 == 3
        cs3 = COSMO.Constraint([1.0 2.0 0 3.0 0 0], [-3.], COSMO.ZeroSet)

        # y2 == 1
        cs4 = COSMO.Constraint([0 0 0 0 1.0 0], [-1.], COSMO.ZeroSet)

        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2; cs3; cs4])

        res = COSMO.optimize!(model)
        @test res.status == :Solved
        @test isapprox(res.obj_val, -1.8458, atol=1e-3)

        # using JuMP
        # m = Model(with_optimizer(COSMO.Optimizer))
        # @variable(m, x1)
        # @variable(m, y)
        # @variable(m, z1)
        # @variable(m, x2)
        # @variable(m, y2 == 1.)
        # @variable(m, z2)

        # @objective(m, Max, z1 + z2)
        # @constraint(m, x1 + 2 * y + 3 * x2 == 3)
        # @constraint(m, [x1; y; z1] in MOI.PowerCone(0.6))
        # @constraint(m, [x2; y2; z2] in MOI.PowerCone(0.1))
        # status = JuMP.optimize!(m)

    end

    @testset "Primal infeasible" begin
        # max  z
        # s.t. x^0.8 * y^0.2 >= z
        #       x = y == 1
        #       z = 2

        P = spzeros(3, 3)
        q = [0; 0; -1.0]

        cs1 = COSMO.Constraint(Matrix(1.0I, 3, 3), zeros(3), COSMO.PowerCone(0.8))
        cs2 = COSMO.Constraint(Matrix(1.0I, 3, 3), [-1.; -1.; -2.], COSMO.ZeroSet)

        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2])

        res = COSMO.optimize!(model)

        @test res.status == :Primal_infeasible
    end

    @testset "Dual infeasible" begin
        # min z
        # s.t. x^0.8 * y^0.2 >= z

        P = spzeros(3, 3)
        q = [0; 0; 1.0]

        cs1 = COSMO.Constraint(Matrix(1.0I, 3, 3), zeros(3), COSMO.PowerCone(0.8))

        model = COSMO.Model()
        assemble!(model, P, q, cs1)

        res = COSMO.optimize!(model)
        @test res.status == :Dual_infeasible
    end

end

@testset "Dual power cone problems" begin
    @testset "Feasible" begin
        # max  z
        # s.t. (x/0.8)^0.8 * (y/0.2)^0.2 >= z
        #       x == 0.8,  y == 0.2
        P = spzeros(3, 3)
        q = [0.; 0.; -1]

        # (x, y, z) in K^*_pow(0.8)
        A1 = sparse(Matrix(1.0I, 3, 3))
        b1 = zeros(3)
        cs1 = COSMO.Constraint(A1, b1, COSMO.DualPowerCone(0.8))

        # x == 0.8, y == 0.2
        cs2 = COSMO.Constraint([1. 0 0; 0 1 0], [-0.8; -0.2], COSMO.ZeroSet)

        model = COSMO.Model()
        assemble!(model, P, q, [cs1; cs2])

        res = COSMO.optimize!(model)
        @test res.status == :Solved
        @test isapprox(res.obj_val, -1., atol=1e-3)

    end
end


nothing
