# Test file to check the projections

workspace()
include("../ossdp.jl")
using OSSDPSolver, Base.Test

rng = MersenneTwister(1234)

nn = 1000

function isNumericallyPosDef(X,eps)
  F = eigfact(X)
  if size(find(x-> x<eps,F[:values]), 1) == 0
    return true
  else
    return false
  end

end
@testset "Test projection onto pos-def cone" begin

  for iii = 1:nn

    # create random matrix
    dim = rand(rng,2:100)
    A = rand(rng,dim,dim)
    A = full(Symmetric(A))
    # project matrix onto positive semidefinite cone
    a = vec(A)
    X = project_sdcone(a,dim)

    # check positive semidefiniteness
    @test isNumericallyPosDef(reshape(X,dim,dim),-1e-10) == true

  end

end