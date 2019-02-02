
using COSMO, Random, Test
rng = Random.MersenneTwister(12345)

include("./UnitTests/COSMOTestUtils.jl")

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
  include("./UnitTests/algebra.jl")
  include("./UnitTests/splitvector.jl")
  include("./UnitTests/interface.jl")
  include("./UnitTests/moi_wrapper.jl")
end
nothing