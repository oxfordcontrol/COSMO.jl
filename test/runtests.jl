include("../../src/QOCS.jl")

using QOCS

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

rng = MersenneTwister(12345)

@testset "All Unit Tests" begin

  include("./UnitTests/simple.jl")
  include("./UnitTests/qp-lasso.jl")
  include("./UnitTests/socp-lasso.jl")
  include("./UnitTests/closestcorr.jl")
  include("./UnitTests/InfeasibilityTests/runTests.jl")

end
