using COSMO, Random, Test, Pkg
rng = Random.MersenneTwister(12345)




#include("../COSMOTestUtils.jl")

@testset "Acceleration test set" begin

  include("./adaptive_rho_acc_restarts.jl")
  include("./max_rho_adaption.jl")
  include("./anderson_accelerator.jl") 
  # include("./delayed_start.jl")
  # include("./safeguarding.jl")

end
nothing
