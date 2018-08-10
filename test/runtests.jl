include("../../src/QOCS.jl")

using QOCS
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

rng = MersenneTwister(12345)

@testset "All tests" begin

  include("./simple.jl")
  include("./qp-lasso.jl")
  include("./socp-lasso.jl")
  include("./closestcorr.jl")
  include("./InfeasibilityTests/runTests.jl")

end


# write your own tests here
@test 1 == 2
