        # Exponential cone constrained problems

using Test, LinearAlgebra, SparseArrays
using COSMO

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

@testset "Exponential cone problems" begin

for T in UnitTestFloats

    @testset "Exponential cone problems (T = $(T))" begin

        @testset "Feasible" begin

            # solve the following exponential cone problem
            # max  x
            # s.t. y * exp(x / y) <= z
            #      y == 1, z == exp(5)

            P = spzeros(T, 3, 3)
            q = T[-1.; 0; 0]

            # y * exp( x / y) <= z
            A1 = spdiagm(0 => ones(T, 3))
            b1 = zeros(T, 3)

            # y == 1 and z == exp(5) = 148.4131591025766
            A2 = SparseMatrixCSC{T, Int}([0 1. 0; 0 0 1])
            b2 = [-one(T); -exp(T(5))]
            cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)
            cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)

            model = COSMO.Model{T}()
            assemble!(model, P, q, [cs1; cs2], settings = COSMO.Settings{T}(eps_abs = T(1e-4), eps_rel = T(1e-4)))

            res = COSMO.optimize!(model)
            @test res.status == :Solved
            @test isapprox(res.obj_val, T(-5), atol=1e-2)
        end

        # Run infeasibility tests only if precision(T) is high enough
        if precision(T) >= precision(Float64)
            @testset "Primal infeasible 1" begin
                # min  x
                # s.t. y * exp( x / y) <= z
                #       y == 1
                #       z == -1

                P = spzeros(T, 3, 3)
                q = T[1.; 0; 0]

                # y * exp( x / y) <= z
                A1 = spdiagm(0 => ones(T, 3))
                b1 = zeros(T, 3)

                # y == 1
                A2 = T[0 -1. 0]
                b2 = T[-1.]

                # z == -1
                A3 = T[0 0 -1.]
                b3 = T[1]

                cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)
                cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)
                cs3 = COSMO.Constraint(A3, b3, COSMO.ZeroSet)

                model = COSMO.Model{T}()
                assemble!(model, P, q, [cs1; cs2; cs3])

                res = COSMO.optimize!(model)
                @test res.status == :Primal_infeasible
            end

            @testset "Primal infeasible 2" begin
                # min  x
                # s.t. y * exp( x / y) <= z - 0.2
                #      y * exp( x / y) >= z + 0.2

                P = spzeros(T, 3, 3)
                q = T[1.; 0; 0]

                # y * exp( x / y) <= z - 0.2
                A1 = spdiagm(0 => ones(T, 3))
                b1 = T[0.; 0; -0.2]

                # y * exp( x / y) >= z + 0.3
                A2 = -spdiagm(0 => ones(T, 3))
                b2 = T[0.; 0; -0.3]

                cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)
                cs2 = COSMO.Constraint(A2, b2, COSMO.ExponentialCone)


                model = COSMO.Model{T}()
                assemble!(model, P, q, [cs1; cs2])

                res = COSMO.optimize!(model)
                @test res.status == :Primal_infeasible
            end

            @testset "Dual infeasible" begin
                # max z
                # s.t. y * exp(x / y) <= z

                P = spzeros(T, 3, 3)
                q = T[0; 0; -1.0]

                # y * exp( x / y) <= z
                A1 = spdiagm(0 => ones(T, 3))
                b1 = zeros(T, 3)

                cs1 = COSMO.Constraint(A1, b1, COSMO.ExponentialCone)

                model = COSMO.Model{T}()
                assemble!(model, P, q, cs1)

                res = COSMO.optimize!(model)
                @test res.status == :Dual_infeasible
            end
        end
    end

    @testset "Dual Exponential cone problems (T = $(T))" begin
        if precision(T) >= precision(Float64)
        @testset "Feasible" begin

            # solve the following dual exponential cone problem
            # min  y
            # s.t. -x * exp(y / x) <= exp(1) * z
            #      x == -1, z == exp(5)

            P = spzeros(T, 3, 3)
            q = T[0; 1.; 0]

            # -x * exp( y / x) <= z
            A1 = spdiagm(0 => ones(T, 3))
            b1 = zeros(T, 3)
            cs1 = COSMO.Constraint(A1, b1, COSMO.DualExponentialCone)

            # y == 1 and z == exp(5)
            A2 = SparseMatrixCSC{T, Int}([1. 0 0; 0 0 1])
            b2 = [one(T); -exp(T(5))]
            cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)

            model = COSMO.Model{T}()
            assemble!(model, P, q, [cs1; cs2], settings = COSMO.Settings{T}())

            res = COSMO.optimize!(model)
            @test res.status == :Solved
            @test isapprox(res.obj_val, T(-6), atol=1e-3)
        end
    end
        # Run infeasibility tests only if precision(T) is high enough
        if precision(T) >= precision(Float64)
            @testset "Primal Infeasible" begin

                # solve the following dual exponential cone problem
                # min  u + v + w
                # s.t. -u * exp(v / u) <= exp(1) * w
                #      u == 1, v == 2

                P = spzeros(T, 3, 3)
                q = ones(T, 3)

                # -u * exp( v / u) <= exp(1) * w
                A1 = spdiagm(0 => ones(T, 3))
                b1 = zeros(T, 3)
                cs1 = COSMO.Constraint(A1, b1, COSMO.DualExponentialCone)

                #  u == 1 and b == 2
                A2 = SparseMatrixCSC{T, Int}([1. 0 0; 0 1 0])
                b2 = T[-1.; -2]
                cs2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet)

                model = COSMO.Model{T}()
                assemble!(model, P, q, [cs1; cs2], settings = COSMO.Settings{T}())

                res = COSMO.optimize!(model)
                @test res.status == :Primal_infeasible

            end
        end
    end
end
end
nothing
