# Test functions related to the reduced clique graph and the merging strategy that uses it
using COSMO, Test, LinearAlgebra, Random


# We use the example graph from the paper in Habib, Stacho - Reduced clique graphs of chordal graphs(2011) - Fig. 1
snd = [Set([4, 5]); Set([1, 4, 6]); Set([1, 7]); Set([1, 8]); Set([1, 3, 4]); Set([1,2,3]); Set([2, 3, 9]); Set([3, 4, 11]); Set([3, 10])]
sep = [Set([1, 3]); Set([1, 4]); Set([2, 3]); Set([3, 4]); Set([1]); Set([3]); Set([4])]

rows, cols = COSMO.compute_reduced_clique_graph!(sep, snd)

edges = Tuple{Int64, Int64}[]
for k = 1:length(rows)
    push!(edges, (rows[k], cols[k]))
end
# compare this to the known reference solution
edges_ref = [(2, 1), (5, 1), (8, 1), (9, 8), (9, 5), (9, 7), (7, 6), (6, 4), (5, 4), (4, 2), (4, 3), (3, 2), (5, 3), (6, 3), (9, 6), (8, 5), (5, 2), (6, 5)]
permissible_ref = falses(length(edges_ref))
permissible_ref[[7, 11, 16, 17, 18]] .= true
@testset "Reduced clique graph" begin

    edges_contained = true
    for edge in edges
        if edge ∉ edges_ref
        edges_contained = false
        break
    end
    end
    @test edges_contained


    # let's make a supernodetree structure to test merging related functions
    N = 11
    t = COSMO.SuperNodeTree(deepcopy(snd), N)
    t.sep = sep
    COSMO.initialise!(t, t.strategy)

    # check whether an edge is permissible
    permissible = falses(length(edges))

    for (k, edge) in enumerate(edges)
        permissible[k] = COSMO.ispermissible(edge, t.strategy.adjacency_table, t.snd)
    end
    permissible_edges_ref = edges_ref[permissible_ref]
    for permissible_edge in edges[permissible]
        @test permissible_edge in permissible_edges_ref
    end


    cand, weight = COSMO.traverse(t, t.strategy)
    # merge the following cliques
    cand = [5, 2];

    @test !COSMO.evaluate(t, t.strategy, cand)
    COSMO.merge_two_cliques!(t, cand, t.strategy)
    @test isempty(t.snd[2])
    @test t.snd[5] == Set([1,3,4,6])

    # update the graph information
    COSMO.update_strategy!(t.strategy, t, cand, true)

    @test !haskey(t.strategy.adjacency_table, 2)
    @test mapreduce(x -> 2 ∈ x, +,  t.strategy.adjacency_table) == 0 #check that all 2's are deleted

    # consider a new merge example: 7 and 6
    N = 11
    t = COSMO.SuperNodeTree(deepcopy(snd), N)
    t.sep = sep
    COSMO.initialise!(t, t.strategy)
    cand = COSMO.traverse(t, t.strategy)
    cand = [7, 6];

    @test !COSMO.evaluate(t, t.strategy, cand)
    COSMO.merge_two_cliques!(t, cand, t.strategy)
    @test isempty(t.snd[6])
    @test t.snd[7] == Set([1,2,3,9])
    COSMO.update_strategy!(t.strategy, t, cand, true)
    @test !haskey(t.strategy.adjacency_table, 6)
    @test mapreduce(x -> 6 ∈ x, +,  t.strategy.adjacency_table) == 0 #check that all 6's are deleted

    # Check the recomputation step
    COSMO.clique_intersections!(t.strategy.edges, t.snd)
    COSMO.kruskal!(t.strategy.edges, t.num)
    COSMO.determine_parent_cliques!(t.snd_par, t.snd_child, t.snd, t.post, t.strategy.edges)
    t.snd_post = COSMO.post_order(t.snd_par, t.snd_child, t.num)
    t.sep = [Set{Int}() for i=1:length(t.snd)]
    COSMO.split_cliques!(t.snd, t.sep, t.snd_par, t.snd_post, t.num)

    @test check_clique_tree(t.snd, t.sep, t.snd_par, t.snd_child)
end #@testset
