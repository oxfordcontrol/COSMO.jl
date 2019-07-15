# Unit tests to check that the implemented merge strategy from CHOMPACK works as expected
using COSMO, SparseArrays, LinearAlgebra, Test


function get_example_tree(strategy)

  # define example tree with known merge order
  snd1 = [15; 16; 17]; sep1 = [];
  snd2 = [5; 9]; sep2 = [15; 16];
  snd3 = [3; 4]; sep3 = [5; 15];
  snd4 = [1]; sep4 = [3];
  snd5 = [2]; sep5 = [3; 4];
  snd6 = [6]; sep6 = [9; 16];
  snd7 = [7; 8]; sep7 = [9; 15];
  snd8 = [12; 13; 14]; sep8 = [16; 17];
  snd9 = [10; 11]; sep9 = [13; 14; 17];

  supernodes = [snd1, snd2, snd3, snd4, snd5, snd6, snd7, snd8, snd9]
  seperators = [sep1, sep2, sep3, sep4, sep5, sep6, sep7, sep8, sep9]
  parents = [0; 1; 2; 3; 3; 2; 2; 1; 8]
  snd_post = collect(9:-1:1)
  post = collect(1:1:17)
  return COSMO.SuperNodeTree(supernodes, parents, snd_post, seperators, strategy, post = post)
end

function get_example_graph(strategy)
  t = get_example_tree(strategy)
  @. t.snd_par = -1
  @. empty!(t.snd_child)

  # add separators to snds
  for (i, snd) in enumerate(t.snd)
    push!(snd, t.sep[i]...)
    sort!(snd)
    t.sep[i] = []
  end
  return t
end


function complexity_savings(t, c1, c2)
  dim_c1 = length(t.snd[c1])
  dim_c2 = length(t.snd[c2])
  dim_u = length(union(t.snd[c1], t.snd[c2]))
  return dim_c1^3 + dim_c2^3 - dim_u^3
end

# Correct solution for adjacency matrix with edge scores
function get_adjacency_matrix(t)
  rows = [2; 3; 6; 7; 8; 9;
      3; 6; 7; 8;
      4; 5; 7;
      5;
      7; 8;
      9]
  cols = [ones(Int64, 6); 2*ones(Int64, 4); 3*ones(Int64, 3); 4; 6*ones(Int64, 2); 8]
  N = length(rows)
  edge_score = zeros(N)
  for iii = 1:N
    edge_score[iii] = complexity_savings(t, rows[iii], cols[iii])
  end
  A = sparse(rows, cols, edge_score, 9, 9)
  return A
end

# correct solution for clique graph with intersection weights
function get_intersection_matrix()
  rows = [2; 3; 6; 7; 8; 9; 3; 6; 7; 8; 4; 5; 7; 5; 7; 8; 9]
  cols = [ones(Int64, 6); 2 * ones(Int64, 4); 3 * ones(Int64, 3); 4; 6*ones(Int64, 2); 8]
  vals = [2.0; 1; 1; 1; 2; 1; 2; 2; 2; 1; 1; 2; 1; 1; 1; 1; 3]
  A = sparse(rows, cols, vals, 9, 9)
end

@testset "Clique merging (CHOMPACK merge strategy)" begin


strategy = COSMO.TreeTraversalMerge()
t = get_example_tree(strategy)
COSMO.initialise!(t, strategy)

# Check first  merge_cand: (1, 2)
cand = COSMO.traverse(t, strategy)
dim_par_snd, dim_par_sep = COSMO.clique_dim(t, cand[1])
dim_clique_snd, dim_clique_sep = COSMO.clique_dim(t, cand[2])
@test cand == [1; 2] && dim_clique_snd == 2 && dim_clique_sep == 2 && dim_par_snd == 3 && dim_par_sep == 0
f_i = COSMO.fill_in(dim_clique_snd, dim_clique_sep, dim_par_snd, dim_par_sep)
s = COSMO.max_snd_size(dim_clique_snd, dim_par_snd)
@test f_i == 2 && s == 3
@test COSMO.evaluate(t, strategy, cand)
ordered_cand = COSMO.merge_two_cliques!(t, cand, strategy)
@test isempty(t.snd[2]) && t.snd[1] == [15; 16; 17; 5; 9]
@test isempty(t.sep[2]) && isempty(t.sep[1])
@test t.num == 8 && t.snd_par[2] == -1 && t.snd_par[3] == t.snd_par[6]  == t.snd_par[7] == 1
@test t.snd_child[1] == [8; 3; 6; 7] && isempty(t.snd_child[2])


# Let's go throught the whole tree and check if the merge log is as expected
strategy = COSMO.TreeTraversalMerge()
t = get_example_tree(strategy)
COSMO.merge_cliques!(t, strategy)
# considered clique pairs
@test t.merge_log.clique_pairs == [1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 8 9]
# decision if merge should be done
@test t.merge_log.decisions == [true; true; true; true; true; false; false; true]
# number of merges
@test t.merge_log.num == 6

end


@testset "Clique merging (Clique graph based merge strategy)" begin


strategy = COSMO.PairwiseMerge(edge_score = COSMO.ComplexityScore())
t = get_example_graph(strategy)
COSMO.initialise!(t, strategy)
# since all values are < 0, this merge strategy wouldn't merge at all
@test strategy.edges == get_adjacency_matrix(t)

t = get_example_graph(strategy)
COSMO.merge_cliques!(t, strategy)
# number of merges
@test t.merge_log.num == 0


# recompute clique tree from clique graph
# compute clique graph edge weights ( |C_i ∩ C_j)
COSMO.clique_intersections!(t.strategy.edges, t.snd)
@test t.strategy.edges == get_intersection_matrix()
# Compute maximum weight spanning tree
COSMO.kruskal!(t.strategy.edges)
@test length(findall(x -> x < 0, t.strategy.edges)) == 8

COSMO.determine_parent_cliques!(t.snd_par, t.snd_child, t.snd, t.post, t.strategy.edges)
# check that tree structure is as expected
@test t.snd_par == [0; 1; 2; 3; 3; 2; 2; 1; 8]
t.snd_post = COSMO.post_order(t.snd_par, t.snd_child)

COSMO.split_cliques!(t.snd, t.sep, t.snd_par, t.snd_post, t.num)


I = [3; 4; 5; 15;3; 4; 4; 5; 15; 5; 15; 9; 15; 16; 9; 16; 8; 9; 15; 9; 15; 15; 16; 11; 13; 14; 17; 13; 14; 17; 13; 14; 16; 17; 14; 16; 17; 16; 17; 16; 17; 17]
J = [ones(Int64, 4); 2 * ones(Int64, 2); 3*ones(Int64, 3); 4*ones(Int64, 2); 5*ones(Int64, 3); 6*ones(Int64, 2); 7*ones(Int64, 3); 8*ones(Int64, 2); 9*ones(Int64, 2); 10*ones(Int64, 4); 11*ones(Int64, 3);
    12*ones(Int64, 4); 13*ones(Int64, 3); 14*ones(Int64, 2); 15*ones(Int64, 2); 16 ]

L = sparse(I, J, ones(length(I)), 17, 17)

strategy = COSMO.PairwiseMerge(edge_score = COSMO.ComplexityScore())
t = COSMO.SuperNodeTree(L, COSMO.TreeTraversalMerge())




end