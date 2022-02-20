# Unit tests to check that the implemented merge strategy from CHOMPACK works as expected
using COSMO, SparseArrays, LinearAlgebra, Test


function get_example_tree(strategy)

  # define example tree with known merge order
  snd1 = Set([15, 16, 17]); sep1 = Set{Int}();
  snd2 = Set([5, 9]); sep2 = Set([15; 16]);
  snd3 = Set([3, 4]); sep3 = Set([5; 15]);
  snd4 = Set([1]); sep4 = Set([3]);
  snd5 = Set([2]); sep5 = Set([3; 4]);
  snd6 = Set([6]); sep6 = Set([9; 16]);
  snd7 = Set([7, 8]); sep7 = Set([9; 15]);
  snd8 = Set([12, 13, 14]); sep8 = Set([16; 17]);
  snd9 = Set([10, 11]); sep9 = Set([13; 14; 17]);

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
    union!(snd, t.sep[i])
    empty!(t.sep[i])
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
  rows = [2;  8; #1
  3; 6; 7; #2
  4; 5; #3
   5; #4
  9; #8
  ]
  cols = [ones(Int64, 2); 2*ones(Int64, 3); 3*ones(Int64, 2);  4; 8]
  N = length(rows)
  edge_weight = zeros(N)
  for iii = 1:N
    edge_weight[iii] = complexity_savings(t, rows[iii], cols[iii])
  end
  A = sparse(rows, cols, edge_weight, 9, 9)
  return A
end

# correct solution for clique graph with intersection weights
function get_intersection_matrix()
  rows = [2; 3; 6; 7; 8; 9; 3; 6; 7; 8; 4; 5; 7; 5; 7; 8; 9]
  cols = [ones(Int64, 6); 2 * ones(Int64, 4); 3 * ones(Int64, 3); 4; 6*ones(Int64, 2); 8]
  vals = [2.0; 1; 1; 1; 2; 1; 2; 2; 2; 1; 1; 2; 1; 1; 1; 1; 3]
  A = sparse(rows, cols, vals, 9, 9)
end
@testset "Clique merging Example" begin

  @testset "Clique merging (Parent-child merge strategy)" begin


    strategy = COSMO.ParentChildMerge()
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
    @test isempty(t.snd[2]) && t.snd[1] == Set([15; 16; 17; 5; 9])
    @test isempty(t.sep[2]) && isempty(t.sep[1])
    @test t.num == 8 && t.snd_par[2] == -1 && t.snd_par[3] == t.snd_par[6]  == t.snd_par[7] == 1
    @test t.snd_child[1] == Set([8; 3; 6; 7]) && isempty(t.snd_child[2])


    # Let's go throught the whole tree and check if the merge log is as expected
    strategy = COSMO.ParentChildMerge()
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

    strategy = COSMO.CliqueGraphMerge(edge_weight = COSMO.ComplexityWeight())
    t = get_example_tree(strategy)
    t.snd = union.(t.snd, t.sep)

    COSMO.initialise!(t, t.strategy)

    # since all values are < 0, this merge strategy wouldn't merge at all
    @test strategy.edges == get_adjacency_matrix(t)

    t = get_example_tree(strategy)
    t.snd = union.(t.snd, t.sep)
    COSMO.merge_cliques!(t, strategy)
    # number of merges
    @test t.merge_log.num == 0
    @test t.snd_par == [0; 1; 2; 3; 3; 2; 2; 1; 8]
  end
end
