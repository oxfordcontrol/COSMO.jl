using COSMO, Test, LinearAlgebra, SparseArrays, Random

@testset "Infeasibility" begin
  include("primal_infeasible_1.jl")
  #include("primal_infeasible_2.jl")
  include("primal_infeasible_3.jl")
  include("dual_infeasible_1.jl")
  include("dual_infeasible_2.jl")
  # include("infeasible_SDPLib.jl")
end
