
using QOCS, Test, SparseArrays, Random, LinearAlgebra


rng = Random.MersenneTwister(1872381)

P_int = 4
q_int = 2

P_Float = 4.
q_Float = 2.

P_UInt = UInt64(4)
q_UInt = UInt64(2)

P_mat = rand(rng,-10:10,10,10)
q_mat = rand(rng,-10:10,10,1)


P_sparse = sprand(rng,10,10,0.4)
q_sparse = rand(rng,1:10,10,1)

A_vec = [1 2 3 4]
b_vec = 1

# A_mat = rand(rng,10,10)
# b_mat = sparse(rand(rng,10,1))
scalar_c = QOCS.Constraint(1,2,QOCS.Zeros())
constr = QOCS.Constraint(Matrix(1.0I,10,10),rand(rng,10),QOCS.Zeros())


@testset "Model" begin
    model = QOCS.Model()
    @test assemble!(model,P_int,q_int,[scalar_c]) == nothing

    model = QOCS.Model()
    @test assemble!(model,P_Float,q_Float,[scalar_c]) == nothing

    model = QOCS.Model()
    @test assemble!(model,P_UInt,q_UInt,[scalar_c]) == nothing

    model = QOCS.Model()
    @test assemble!(model,P_mat,q_mat,[constr]) == nothing

    model = QOCS.Model()
    @test assemble!(model,P_sparse,q_sparse,[constr]) == nothing

end
nothing



