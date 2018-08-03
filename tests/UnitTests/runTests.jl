include("../../src/Solver.jl")

using OSSDP, OSSDPTypes, Base.Test

rng = MersenneTwister(12345)

@testset "All tests" begin

  include("./simple.jl")
  include("./qp-lasso.jl")
  include("./socp-lasso.jl")
  include("./closestcorr.jl")
  include("./InfeasibilityTests/runTests.jl")

end
