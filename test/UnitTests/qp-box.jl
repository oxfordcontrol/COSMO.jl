# Box constrained problems

using Test, LinearAlgebra, SparseArrays, Random
using COSMO

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end
@testset "Box Problems" begin

for T in UnitTestFloats

    @testset "Box Problems (T = $(T))" begin
        @testset "Box Problems - feasible" begin
            A = SparseMatrixCSC{T, Int64}([1. 0; 0. 1.])
            b = T[0.;0.]

            P = Matrix(Diagonal(ones(T, 2)))
            q = T[1.;-1.]
            l = T[0.;0]
            u = T[1.;1.]
            constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
            constraints = [constraint]

            model = COSMO.Model{T}()
            assemble!(model, P, q, constraints, settings = COSMO.Settings{T}(alpha = 1.4))

            res = COSMO.optimize!(model)
            @test res.status == :Solved
            @test isapprox(res.obj_val, -0.5, atol=1e-5)
        end
        nothing

        @testset "Box Problems - primal infeasible" begin
            A = SparseMatrixCSC{T, Int64}([1. 0; 1. 0.])
            b = T[2.;0.]

            P = Matrix(Diagonal(ones(T, 2)))
            q = T[1.;-1.]
            l = T[0.;0]
            u = T[1.;1.]
            constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
            constraints = [constraint]

            model = COSMO.Model{T}()
            assemble!(model, P, q, constraints)

            res = COSMO.optimize!(model)
            @test res.status == :Primal_infeasible
        end

        @testset "Box Problems - primal infeasible" begin
            A = SparseMatrixCSC{T, Int64}([1. 0; 1. 0.])
            b = T[0.;0.]

            P = Matrix(Diagonal(ones(T, 2)))
            q = T[1.;-1.]
            l = T[0.;2.]
            u = T[1.;3.]
            constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
            constraints = [constraint]

            model = COSMO.Model{T}()
            assemble!(model, P, q, constraints)

            res = COSMO.optimize!(model)
            @test res.status == :Primal_infeasible
        end

        @testset "Box Problems - dual infeasible (unscaled)" begin
            A = SparseMatrixCSC{T, Int64}([1. 0; 0. 1.])
            b = T[1.;1.]

            P = Matrix(Diagonal(zeros(T, 2)))
            q = T[1.;1.]
            l = T[0.;-Inf]
            u = T[1.;3.]
            constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
            constraints = [constraint]

            model = COSMO.Model{T}()
            settings = COSMO.Settings{T}(check_infeasibility=20,scaling=0)
            assemble!(model, P, q, constraints,settings  = settings)

            res = COSMO.optimize!(model)
            @test res.status == :Dual_infeasible
        end

        @testset "Box Problems - dual infeasible (scaled)" begin
            A = SparseMatrixCSC{T, Int64}([1. 0; 0. 1.])
            b = T[1.;1.]

            P = Matrix(Diagonal(zeros(T, 2)))
            q = T[1.;1.]
            l = T[0.;-Inf]
            u = T[1.;3.]
            constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
            constraints = [constraint]
            model = COSMO.Model{T}()
            settings = COSMO.Settings{T}(check_infeasibility=40,scaling=10)
            assemble!(model, P, q, constraints,settings = settings)

            res = COSMO.optimize!(model)
            @test res.status == :Dual_infeasible
        end
    end
end
end
nothing
