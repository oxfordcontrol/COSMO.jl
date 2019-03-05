# Box constrained problems

using Test, LinearAlgebra, SparseArrays, Random
using COSMO


@testset "Box Problems - feasible" begin
    A = sparse([1. 0; 0. 1.])
    b = [0.;0.]

    P = Matrix(Diagonal(ones(2)))
    q = [1.;-1.]
    l = [0.;0]
    u = [1.;1.]
    constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
    constraints = [constraint]

    model = COSMO.Model()
    assemble!(model, P, q, constraints)

    res = COSMO.optimize!(model)
    @test res.status == :Solved
    @test isapprox(res.obj_val, -0.5, atol=1e-5)
end
nothing

@testset "Box Problems - primal infeasible" begin
    A = sparse([1. 0; 1. 0.])
    b = [2.;0.]

    P = Matrix(Diagonal(ones(2)))
    q = [1.;-1.]
    l = [0.;0]
    u = [1.;1.]
    constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
    constraints = [constraint]

    model = COSMO.Model()
    assemble!(model, P, q, constraints)

    res = COSMO.optimize!(model)
    @test res.status == :Primal_infeasible
end

@testset "Box Problems - primal infeasible" begin
    A = sparse([1. 0; 1. 0.])
    b = [0.;0.]

    P = Matrix(Diagonal(ones(2)))
    q = [1.;-1.]
    l = [0.;2.]
    u = [1.;3.]
    constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
    constraints = [constraint]

    model = COSMO.Model()
    assemble!(model, P, q, constraints)

    res = COSMO.optimize!(model)
    @test res.status == :Primal_infeasible
end

@testset "Box Problems - dual infeasible (unscaled)" begin
    A = sparse([1. 0; 0. 1.])
    b = [1.;1.]

    P = Matrix(Diagonal(zeros(2)))
    q = [1.;1.]
    l = [0.;-Inf]
    u = [1.;3.]
    constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
    constraints = [constraint]

    model = COSMO.Model()
    settings = COSMO.Settings(check_infeasibility=20,scaling=0)
    assemble!(model, P, q, constraints,settings  = settings)

    res = COSMO.optimize!(model)
    @test res.status == :Dual_infeasible
end

@testset "Box Problems - dual infeasible (scaled)" begin
    A = sparse([1. 0; 0. 1.])
    b = [1.;1.]

    P = Matrix(Diagonal(zeros(2)))
    q = [1.;1.]
    l = [0.;-Inf]
    u = [1.;3.]
    constraint = COSMO.Constraint(A,b, COSMO.Box(l,u))
    constraints = [constraint]

    model = COSMO.Model()
    settings = COSMO.Settings(check_infeasibility=40,scaling=10)  
    assemble!(model, P, q, constraints,settings = settings)

    res = COSMO.optimize!(model)
    @test res.status == :Dual_infeasible
end



nothing
