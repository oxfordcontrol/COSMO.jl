using COSMO, Random, Test

include("./UnitTests/COSMOTestUtils.jl")

@testset "All Unit Tests" begin

  # Split tests into native interface unit tests and MOI unit tests
  include("./run_cosmo_tests.jl")
  # include("./UnitTests/moi_wrapper.jl")

end
nothing
