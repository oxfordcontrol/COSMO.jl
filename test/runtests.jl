
using QOCS, Random

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

rng = Random.MersenneTwister(12345)

@testset "All Unit Tests" begin

  include("./UnitTests/simple.jl")
  include("./UnitTests/sets.jl")
  include("./UnitTests/constraints.jl")
  include("./UnitTests/model.jl")
  include("./UnitTests/qp-lasso.jl")
  include("./UnitTests/socp-lasso.jl")
  include("./UnitTests/closestcorr.jl")
  include("./UnitTests/print.jl")
  include("./UnitTests/InfeasibilityTests/runTests.jl")
end
nothing