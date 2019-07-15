# Unit tests for clique merging
using COSMO, SparseArrays, LinearAlgebra, Test


# mutable struct SnodeTree
#   res::Array{Array{Int64, 1},1}
#   sep::Array{Array{Int64, 1},1}
#   par::Array{Int64}
#   post::Array{Int64}
#   child::Array{Array{Int64, 1},1}
#   function SnodeTree(res, sep, par, post)
#     # compute children structure
#     child = COSMO.child_from_par(par)
#     new(res, sep, par, post, child)
#   end
# end

# let's consider the following graph with cliques (seperators | residuals)
#       ---------
#       |  | 1,2 |
#       ---------
#      /         \
# ---------     ---------
# | 1 | 3  |    | 2 | 4  |
# ---------     ---------
res1 = [[1; 2], [3], [4]];
sep1 = [[], [1], [2]];
par1 = [0; 1; 1];
post1 = [1; 2; 3]
ref_edges = spzeros(3, 3)
ref_edges[1, 2] = 1 / 3;
ref_edges[1, 3] = 1 / 3;
ref_edges[3, 2] = 0;


t = SnodeTree(res1, sep1, par1, post1)

@testset "Clique merging (simple example)" begin

  # compute and verify the initial (approximate) "clique graph" layout and edge values (only consider parent-child and sibling-sibling)
  edges = COSMO.initial_clique_graph(t)
  @test edges == ref_edges

  # merge first pair
  cand = COSMO.find_merge_candidates(edges)
  @test cand[1] == 1 && cand[2] == 2

  COSMO.merge_cliques!(t, edges, cand)
  @test t.res == [[1; 2; 3],[], [4]]
  @test t.sep == [[], [], [2]]
  @test t.par == [0, -1, 1]
  COSMO.update!(t, edges, cand)
  ref_edges[1, 2] = 0
  ref_edges[1, 3] = 1 / 4;
  @test edges == ref_edges

  # merge remaining pair
  cand = COSMO.find_merge_candidates(edges)
  COSMO.merge_cliques!(t, edges, cand)
  COSMO.update!(t, edges, cand)
  @test t.res == [[1; 2; 3; 4],[], []]
  @test t.sep == [[], [], []]
  @test t.par == [0, -1, -1]
  ref_edges[1, 3] = 0
  @test edges == ref_edges
end

# lets consider the following graph with cliques (seperators | residuals)
#       ----------------------
#       |  | 1,2,3,6,7,8,9,10 |
#       ----------------------
#      /         \
# ---------     -----------
# | 1 | 5  |    | 1,2 | 4 |
# ---------     -----------
res2 = [[1; 2; 3; 6; 7; 8; 9; 10], [5], [4]];
sep2 = [[], [1], [1;2]];
par2 = [0; 1; 1];
post2 = [1; 2; 3]

t = SnodeTree(res2, sep2, par2, post2)


@testset "Clique merging (sibling merge)" begin


  edges = COSMO.initial_clique_graph(t)

  # merge first pair (siblings)
  cand = COSMO.find_merge_candidates(edges)
  @test cand[1] == 2 && cand[2] == 3

  COSMO.merge_cliques!(t, edges, cand)
  @test t.res == [[1; 2; 3; 6; 7; 8; 9; 10], [], [4; 5]];
  @test t.sep == [[], [], [1; 2]];
  @test t.par == [0, -1, 1]
  COSMO.update!(t, edges, cand)
  @test edges[1, 3] == 0.2
end

# lets consider the following more complex clique tree
# from Vandenberghe: Chordal Graphs and Semidefinite Optimization p. 287 (Fig 4.7, left)


res3 = [[12; 13; 14; 16; 17], [15], [5; 9], [1; 3; 4], [2], [6], [7; 8], [10; 11]];
sep3 = [[], [16; 17], [15; 16], [5; 15], [3; 4], [9; 16], [9; 15], [13; 14; 17]];
par3 = [0; 1; 2; 3; 4; 3; 3; 1];
post3 = collect(1:8)

t = SnodeTree(res3, sep3, par3, post3)


@testset "Clique merging (large example)" begin


  edges = COSMO.initial_clique_graph(t)

  cand = COSMO.find_merge_candidates(edges)
  @test cand[1] == 1 && cand[2] == 8

  COSMO.merge_cliques!(t, edges, cand)
  COSMO.update!(t, edges, cand)
  cand = COSMO.find_merge_candidates(edges)
  @test cand[1] == 2 && cand[2] == 3

  COSMO.merge_cliques!(t, edges, cand)
  COSMO.update!(t, edges, cand)
  cand = COSMO.find_merge_candidates(edges)
  @test cand[1] == 4 && cand[2] == 5


  COSMO.merge_cliques!(t, edges, cand)
  COSMO.update!(t, edges, cand)

  # merge until one node remains
  for i = 1:4
    cand = COSMO.find_merge_candidates(edges)
    COSMO.merge_cliques!(t, edges, cand)
    COSMO.update!(t, edges, cand)
  end

end
