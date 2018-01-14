# Test file to check the projections

workspace()
include("../Projections.jl")
include("../Helper.jl")
using Base.Test, Helper

rng = MersenneTwister(1234)

nn = 1000

@testset "Test projection onto pos-def cone" begin

  for iii = 1:nn

    # create random matrix
    dim = rand(rng,1:100)
    A = rand(rng,dim,dim)
    A = full(Symmetric(A))
    # project matrix onto positive semidefinite cone
    a = vec(A)
    X = Projections.sdcone(a,dim)

    # check positive semidefiniteness
    @test isNumericallyPosSemDef(reshape(X,dim,dim),-1e-10) == true

  end

end