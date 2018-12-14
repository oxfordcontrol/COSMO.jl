
using COSMO, Test, SparseArrays, Random, LinearAlgebra


rng = Random.MersenneTwister(1872381)

P_int = 4
q_int = 2

P_Float = 4.
q_Float = 2.

P_UInt = UInt64(4)
q_UInt = UInt64(2)

P_mat = rand(rng, -10.:10., 10, 10) #float64 types
q_mat = rand(rng, -10.:10., 10, 1)


P_sparse = sprand(rng, 10, 10, 0.4) #float64 types
q_sparse = rand(rng, 1.:10.10, 10, 1)

A_vec = [1 2 3 4]
b_vec = 1


scalar_c = COSMO.Constraint(1., 2., COSMO.ZeroSet)
constr = COSMO.Constraint(Matrix(1.0I, 10, 10), rand(rng, 10), COSMO.ZeroSet)


@testset "Model" begin
     model = COSMO.Model()
     @test assemble!(model,P_int,q_int,[scalar_c]) == nothing
    #
     model = COSMO.Model{typeof(P_Float)}()
    @test assemble!(model,P_Float,q_Float,[scalar_c]) == nothing

    model = COSMO.Model()       #default type
    @test assemble!(model, P_Float, q_Float, [scalar_c]) == nothing

    # model = COSMO.Model{typeof(P_UInt)}()
    # @test assemble!(model,P_UInt,q_UInt,[scalar_c]) == nothing

    model = COSMO.Model()
    @test assemble!(model, P_mat, q_mat, [constr]) == nothing

    model = COSMO.Model()
    @test assemble!(model, P_sparse, q_sparse, [constr]) == nothing

end
nothing
