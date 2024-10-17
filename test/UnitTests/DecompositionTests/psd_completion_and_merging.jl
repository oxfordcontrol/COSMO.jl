# This unit tests checks multiple thigs:
# - different merge strategies should lead to the same results
# - the PSD completion should produce positive semidefinite dual variables Y
# - this also implicetely tests the recomputation of a clique tree from a clique graph when
# a graph based merge strategy is used

using COSMO, Random, Test, LinearAlgebra, SparseArrays
#include("./../COSMOTestUtils.jl")
rng = Random.MersenneTwister(375)



# PROBLEM 1:
#
# An example problem where you can decompose a 5x5 matrix into two
# 4x4 blocks that overlap so much, that you should merge them back into a 5x5 block.
# S = [x x 0 x x
#      x x x x x
#      0 x x x x
#      x x x x x
#.     x x x x x ]
pattern = ones(5, 5)
pattern[1, 3] = 0
pattern = Matrix(Symmetric(pattern, :U))

# PROBLEM 2:
#
# An example with four cliques C1: {1,2,3,4}, C2:{2,3,4,5}, C3: {4,5,6,7,8}, C4:{5,6,7,8,9} (without AMD reordering)
# where |C1 ∩ C2 | = 3 and |C3 ∩ C4| = 4, so idealy you want two merges to happen
# S = [x x x x 0 0 0 0 0
#        x x x x 0 0 0 0
#          x x x 0 0 0 0
#            x x x x x 0
#.             x x x x x
#                x x x x
#                  x x x
#                    x x
#                      x]
pattern2 = ones(9, 9)
pattern2[1, 5:9] .= 0
pattern2[2, 6:9] .= 0
pattern2[3, 6:9] .= 0
pattern2[4, 9] = 0
pattern2 = Matrix(Symmetric(pattern2, :U))


# PROBLEM 3:
# We consider the sparsity pattern from Vandenberghe - Chordal Graphs and semidefinite optimisation (2015) p.279
R = [3; 4; 5; 15;3; 4; 4; 5; 15; 5; 15; 9; 15; 16; 9; 16; 8; 9; 15; 9; 15; 15; 16; 11; 13; 14; 17; 13; 14; 17; 13; 14; 16; 17; 14; 16; 17; 16; 17; 16; 17; 17]
C = [ones(Int64, 4); 2 * ones(Int64, 2); 3*ones(Int64, 3); 4*ones(Int64, 2); 5*ones(Int64, 3); 6*ones(Int64, 2); 7*ones(Int64, 3); 8*ones(Int64, 2); 9*ones(Int64, 2); 10*ones(Int64, 4); 11*ones(Int64, 3);
    12*ones(Int64, 4); 13*ones(Int64, 3); 14*ones(Int64, 2); 15*ones(Int64, 2); 16 ]
pattern3 = sparse(R, C, ones(length(C)), 17, 17)  + Matrix(1.0I, 17, 17)
pattern3 = Symmetric(pattern3, :L)


problem_patterns = [pattern, pattern2, pattern3]


@testset "PSD Completion with different merge strategies" begin
  @testset "Pattern Nr $(k)" for k = 1:length(problem_patterns)

  P, q, At, bt, Ct = feasible_sdp_with_pattern(rng, problem_patterns[k])


  # Run four test cases (no decomposition and decomposition with three different merge strategies)
  cases = [
    COSMO.Settings(decompose = false);
    COSMO.Settings(decompose = true, complete_dual = true, merge_strategy = COSMO.NoMerge);
    COSMO.Settings(decompose = true, complete_dual = true, merge_strategy = COSMO.ParentChildMerge);
    COSMO.Settings(decompose = true, complete_dual = true, merge_strategy = COSMO.CliqueGraphMerge);
    ]
  results = Array{COSMO.Result{Float64}}(undef, 4);

  for i = 1:4
    model = COSMO.Model()
    settings = cases[i]
    COSMO.set!(model, P, q, At, bt, Ct, settings)
    res = COSMO.optimize!(model);
    results[i] = res
  end

    @test abs(results[1].obj_val - results[2].obj_val) < 1e-3 && abs(results[1].obj_val - results[3].obj_val) < 1e-3 && abs(results[1].obj_val - results[4].obj_val) < 1e-3
    @test norm(results[1].x - results[2].x, Inf) < 1e-3 && norm(results[1].x - results[3].x, Inf) < 1e-3 && norm(results[1].x - results[4].x, Inf) < 1e-3
    # check that completed matrices are positive semidefinite
    Y1 = matrixify(results[1].y)
    Y2 = matrixify(results[1].y)
    Y3 = matrixify(results[1].y)
    Y4 = matrixify(results[1].y)

    @test minimum(eigen(Y1).values) > -1e-6
    @test minimum(eigen(Y2).values) > -1e-6
    @test minimum(eigen(Y3).values) > -1e-6
    @test minimum(eigen(Y4).values) > -1e-6
  end

end
