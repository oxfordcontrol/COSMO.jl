using COSMO, Random, Test, Pkg
rng = Random.MersenneTwister(12345)




include("./UnitTests/COSMOTestUtils.jl")

@testset "COSMO native tests" begin

  include("./UnitTests/simple.jl")
  include("./UnitTests/sets.jl")
  include("./UnitTests/constraints.jl")
  include("./UnitTests/kktsolver.jl")
  include("./UnitTests/model.jl")
  include("./UnitTests/qp-box.jl")
  include("./UnitTests/socp-lasso.jl")
  include("./UnitTests/closestcorr.jl")
  include("./UnitTests/exp_cone.jl")
  include("./UnitTests/pow_cone.jl")
  include("./UnitTests/print.jl")
  include("./UnitTests/InfeasibilityTests/runTests.jl")
  include("./UnitTests/algebra.jl")
  include("./UnitTests/splitvector.jl")
  include("./UnitTests/interface.jl")
  include("./UnitTests/chordal_decomposition_triangle.jl")
  include("./UnitTests/nuclear_norm_minimization.jl")
  include("./UnitTests/psd_completion.jl")
  include("./UnitTests/psd_completion_and_merging.jl")
  include("./UnitTests/clique_merging_example.jl")

    # optional unittests
  if pkg_installed("Pardiso", "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2")
    include("./UnitTests/options_factory.jl")
  end

end
