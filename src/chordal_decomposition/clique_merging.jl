"""
    AbstractEdgeWeight

A supertype for the weights that are put on the edges of the clique graph.
"""
abstract type AbstractEdgeWeight end

"""
    AbstractComplexityWeight

A supertype for weights that are based on the complexity of clique-related operations in the
solver algorithm.

All subtypes have to define the function:
`compute_complexity_savings(n_1::Int, n_2::Int, n_m::Int, edge_weight::Subtype) --> Float64`

with
  - n_1: `|C_1|`
  - n_2: `|C_2|`
  - n_3: `|C_1 ∪ C_2|`
"""
abstract type AbstractComplexityWeight <: AbstractEdgeWeight end


"""
    ComplexityWeight

An edge weight that approximates the computational savings from merging two cliques
`Ci` and `Cj` by `e(Ci, Cj) = |Ci|^3 + |Cj|^3 - |Ci ∪ Cj|^3`.
"""
struct ComplexityWeight <: AbstractComplexityWeight end


"""
    NoMerge <: AbstractMergeStrategy

A strategy that does not merge cliques.
"""
struct NoMerge <: AbstractTreeBasedMerge end


"""
    CliqueGraphMerge(edge_weight::AbstractEdgeWeight = ComplexityWeight()) <: AbstractGraphBasedMerge

The (default) merge strategy based on the *reduced* clique graph ``\\mathcal{G}(\\mathcal{B}, \\xi)``, for a set of cliques ``\\mathcal{B} = \\{ \\mathcal{C}_1, \\dots, \\mathcal{C}_p\\}`` where the edge set ``\\xi`` is obtained by taking the edges of the union of clique trees.

Moreover, given an edge weighting function ``e(\\mathcal{C}_i,\\mathcal{C}_j) = w_{ij}``, we compute a weight for each edge that quantifies the computational savings of merging the two cliques.
After the initial weights are computed, we merge cliques in a loop:

**while** clique graph contains positive weights:
- select two permissible cliques with the highest weight ``w_{ij}``
- merge cliques ``\\rightarrow`` update clique graph
- recompute weights for updated clique graph

Custom edge weighting functions can be used by defining your own `CustomEdgeWeight <: AbstractEdgeWeight` and a corresponding `edge_metric` method. By default, the `ComplexityWeight <: AbstractEdgeWeight` is used which computes the weight based
on the cardinalities of the cliques: ``e(\\mathcal{C}_i,\\mathcal{C}_j)  = |\\mathcal{C}_i|^3 + |\\mathcal{C}_j|^3 - |\\mathcal{C}_i \\cup \\mathcal{C}_j|^3``.

See also: *Garstka, Cannon, Goulart - A clique graph based merging strategy for decomposable SDPs (2019)*
"""
mutable struct CliqueGraphMerge <: AbstractGraphBasedMerge
  stop::Bool                                  # a flag to indicate that merging should be stopped
  edges::SparseMatrixCSC{Float64, Int}      # the edges and weights of the reduced clique graph
  p::Array{Int, 1}                          # as a workspace variable to store the sorting of weights
  adjacency_table::Dict{Int, Set{Int}}    # a double structure of edges, to allow fast lookup of neighbors
  edge_weight::AbstractEdgeWeight             # used to dispatch onto the correct scoring function
  clique_tree_recomputed::Bool                # a flag to indicate whether a final clique tree has been recomputed from the clique graph
  function CliqueGraphMerge(; edge_weight = ComplexityWeight())
    new(false, spzeros(0, 0),  Int[], Dict{Int, Set{Int}}(), edge_weight, false)
  end
end


"""
    ParentChildMerge(t_fill = 8, t_size = 8) <: AbstractTreeBasedMerge

The merge strategy suggested in *Sun and Andersen - Decomposition in conic optimization with partially separable structure (2014)*.
The initial clique tree is traversed in topological order and a clique ``\\mathcal{C}_\\ell``  is greedily merged to its parent clique
``\\mathcal{C}_{par(\\ell)}`` if at least one of the two conditions are met

  - ``(| \\mathcal{C}_{par(\\ell)}| -| \\eta_\\ell|) (|\\mathcal{C}_\\ell| - |\\eta_\\ell|) \\leq t_{\\text{fill}}`` (fill-in condition)
  - ``\\max \\left\\{ |\\nu_{\\ell}|,  |\\nu_{par(\\ell)}|  \\right\\} \\leq t_{\\text{size}}`` (supernode size condition)
"""
mutable struct ParentChildMerge <: AbstractTreeBasedMerge
  stop::Bool
  clique_ind::Int
  t_fill::Int
  t_size::Int

  function ParentChildMerge(; t_fill = 8, t_size = 8)
    new(false, 2, t_fill, t_size)
  end
end

"Free unneccessary memory after merging completed."
function free_clique_graph!(strategy::AbstractGraphBasedMerge)
  strategy.edges = spzeros(0, 0)
  strategy.p = Int[]
  strategy.adjacency_table = Dict{Int, Set{Int}}()
  return nothing
end

# Main clique merging routine:
# 1. initialise!() - tree and strategy
# 2. traverse() - find the next merge candidates
# 3. evalute() - evaluate wether to merge them or not
# 4. log() - log the decision
# 5. update_strategy!() - update the tree/graph and strategy
function _merge_cliques!(t::SuperNodeTree, strategy::AbstractMergeStrategy)
  initialise!(t, strategy)

  while !strategy.stop
    # find merge candidates

    cand = traverse(t, strategy)
    # evaluate wether to merge the candidates
    do_merge = evaluate(t, strategy, cand)
    if do_merge
      merge_two_cliques!(t, cand, strategy)
    end
    log_merge!(t, do_merge, cand)

    # update strategy information after the merge
    update_strategy!(strategy, t, cand, do_merge)
    t.num == 1 && break
    strategy.stop && break

  end
  return nothing
end

"Given a supernode elimination tree `t` and a merge strategy `strategy`, merge the cliques in the tree."
merge_cliques!(t::SuperNodeTree) = merge_cliques!(t, t.strategy)
merge_cliques!(t::SuperNodeTree, strategy::NoMerge) = nothing

function merge_cliques!(t::SuperNodeTree, strategy::AbstractTreeBasedMerge)
  # call main merge routine
  _merge_cliques!(t, strategy)
  # do some merge strategy specific post-processing of the tree
  Nc = t.num
  # the merging creates empty supernodes and seperators, recalculate a post order for the supernodes
  snd_post = post_order(t.snd_par, t.snd_child, Nc)
  t.snd_post = snd_post

  return nothing
end

function merge_cliques!(t::SuperNodeTree, strategy::AbstractGraphBasedMerge)
  # call main merge routine
  _merge_cliques!(t, strategy)
  # since for now we have a graph, not a tree, a post ordering or a parent structure does not make sense. Therefore just number
  # the non-empty supernodes in t.snd
  t.snd_post = findall(x -> !isempty(x), t.snd)
  t.snd_par = -ones(Int, length(t.snd))

  # recompute a clique tree from the clique graph
  t.num > 1 && clique_tree_from_graph!(t)

  #turn the unordered sets for t.snd and t.sep into orderable arrays
  t.snd = sort.(collect.(t.snd))
  t.sep = sort.(collect.(t.sep))

  # free up memory
  free_clique_graph!(t.strategy)
  return nothing
end

# calculate block sizes (notice: the block sizes are stored in post order)
function calculate_block_dimensions!(t::SuperNodeTree)
  Nc = t.num
  t.nBlk = zeros(Nc)
  for iii = 1:Nc
    c = t.snd_post[iii]
    t.nBlk[iii] = length(t.sep[c]) + length(t.snd[c])
  end
end

" Merge two cliques `cand` of the tree `t` that are in a parent - child relationship."
function merge_child!(t::SuperNodeTree, cand::Array{Int, 1})
  # determine which clique is the parent
  p, ch = determine_parent(t, cand[1], cand[2])

  # merge child's vertex sets into parent's vertex set
  union!(t.snd[p], t.snd[ch])
  empty!(t.snd[ch])
  empty!(t.sep[ch])

  # update parent structure
  for grandch in t.snd_child[ch]
    t.snd_par[grandch] = p
  end
  t.snd_par[ch] = -1 #-1 instead of NaN, effectively remove that entry from the parent list

  # update children structure
  delete!(t.snd_child[p], ch)
  union!(t.snd_child[p], t.snd_child[ch])
  empty!(t.snd_child[ch])

  # decrement number of cliques in tree
  t.num -= 1
  return nothing
end

"Given the clique graph `t` merge the two cliques with indices in `cand`."
function merge_two_cliques!(t::SuperNodeTree, cand::Array{Int, 1}, strategy::AbstractGraphBasedMerge)
  c_1 = cand[1]
  c_2 = cand[2]
  # merge clique c_2 into c_1
  union!(t.snd[c_1], t.snd[c_2])
  empty!(t.snd[c_2])

  # decrement number of cliques in tree
  t.num -= 1

  return nothing
end

"Given the clique tree `t` merge the two cliques with indices in `cand` as parent-child."
merge_two_cliques!(t::SuperNodeTree, cand::Array{Int, 1}, strategy::ParentChildMerge) = merge_child!(t, cand)

# compute the edges and weights of the initial reduced clique graph
function initialise!(t::SuperNodeTree, strategy::CliqueGraphMerge)
  # save("rs1555_before_merging.jld2", "snd", t.snd, "sep", t.sep)
  # compute the edges and intersections of cliques in the reduced clique graph
  rows, cols = compute_reduced_clique_graph!(t.sep, t.snd)

  weights = compute_weights!(rows, cols, t.snd, strategy.edge_weight)

  strategy.edges = sparse(rows, cols, weights, t.num, t.num)
  strategy.p = zeros(Int, length(strategy.edges.nzval))
  strategy.adjacency_table = compute_adjacency_table(strategy.edges, t.num)
  return nothing
end

function initialise!(t::SuperNodeTree, strategy::ParentChildMerge)
  # start with node that has second highest order
  strategy.clique_ind = length(t.snd) - 1
end


"Compute the edge weight between all cliques specified by the edges (rows, cols)."
function compute_weights!(rows::Array{Int, 1}, cols::Array{Int, 1}, snd::Array{Set{Int}, 1}, edge_weight::AbstractComplexityWeight)
  weights = zeros(Float64, length(rows))
  for k = 1:length(rows)
    c_1 = snd[rows[k]]
    c_2 = snd[cols[k]]
    weights[k] = edge_metric(c_1, c_2, edge_weight)
  end
  return weights
end

"Find the next two cliques in the clique graph `t` to merge."
function traverse(t::SuperNodeTree, strategy::CliqueGraphMerge)

   p = strategy.p
   # find edge with highest weight, if permissible return cliques
   edge = max_elem(strategy.edges)
   ispermissible(edge, strategy.adjacency_table, t.snd) && return [edge[1]; edge[2]]
   # else: sort the weights in edges.nzval to find the permutation p
   sortperm!(view(p, 1:length(strategy.edges.nzval)), strategy.edges.nzval, alg = QuickSort, rev = true)

   # try edges with decreasing weight and check if the edge is permissible
  for k = 2:length(strategy.edges.nzval)
    edge = edge_from_index(strategy.edges, p[k])
    if ispermissible(edge, strategy.adjacency_table, t.snd)
      return [edge[1]; edge[2]]
    end
  end

end

" Traverse tree `t` in descending topological order and return parent and clique (root has highest order)."
function traverse(t::SuperNodeTree, strategy::ParentChildMerge)
  c = t.snd_post[strategy.clique_ind]
  return [t.snd_par[c]; c]
end

"Decide whether to merge the two cliques with clique indices `cand`."
function evaluate(t::SuperNodeTree, strategy::ParentChildMerge, cand::Array{Int, 1})
  strategy.stop && return false
  par = cand[1]
  c = cand[2]
  dim_par_snd, dim_par_sep = COSMO.clique_dim(t, par)
  dim_clique_snd, dim_clique_sep = COSMO.clique_dim(t, c)
  return fill_in(dim_clique_snd, dim_clique_sep, dim_par_snd, dim_par_sep) <= strategy.t_fill || max_snd_size(dim_clique_snd, dim_par_snd) <= strategy.t_size
end


function evaluate(t::SuperNodeTree, strategy::CliqueGraphMerge, cand::Array{Int, 1})
  do_merge = (strategy.edges[cand[1], cand[2]] >= 0)

  if !do_merge
    strategy.stop = true
  end
  return do_merge
end

" After a merge attempt, update the strategy information."
function update_strategy!(strategy::ParentChildMerge, t::SuperNodeTree, cand::Array{Int, 1}, do_merge::Bool)
  # try to merge last node of order 1, then stop
  if strategy.clique_ind == 1
    strategy.stop = true
    # otherwise decrement node index
  else
    strategy.clique_ind -= 1
  end
end

"After a merge happened, update the reduced clique graph."
function update_strategy!(strategy::CliqueGraphMerge, t::SuperNodeTree, cand::Array{Int, 1}, do_merge::Bool)

  # After a merge operation update the information of the strategy
  if do_merge

    c_1_ind = cand[1]
    c_removed = cand[2]
    edges = strategy.edges
    n = size(edges, 2)
    adjacency_table = strategy.adjacency_table

    c_1 = t.snd[c_1_ind]
    neighbors = adjacency_table[c_1_ind]
    # neighbors exclusive to the removed clique (and not c1)
    new_neighbors = setdiff(adjacency_table[c_removed], neighbors, c_1_ind)


    # recalculate edge values of all of c_1's neighbors
    for n_ind in neighbors
        if n_ind != c_removed
          neighbor = t.snd[n_ind]
          edges[max(c_1_ind, n_ind), min(c_1_ind, n_ind)] = COSMO.edge_metric(c_1, neighbor, strategy.edge_weight)
        end
    end

    # point edges exclusive to removed clique to "surviving" clique 1
    for n_ind in new_neighbors
        neighbor = t.snd[n_ind]
        edges[max(c_1_ind, n_ind), min(c_1_ind, n_ind)]  = COSMO.edge_metric(c_1, neighbor, strategy.edge_weight)
    end

    # overwrite the weight to any "deleted" edges that still contain a link to c_removed
    strategy.edges[c_removed+1:n, c_removed] .= 0
    strategy.edges[c_removed, 1:c_removed] .= 0
    dropzeros!(edges)

    # update adjacency table in a similar manner
    union!(adjacency_table[c_1_ind], new_neighbors)
    for new_neighbor in new_neighbors
      push!(adjacency_table[new_neighbor], c_1_ind)
    end
    delete!(adjacency_table, c_removed)
    for (key, set) in adjacency_table
      delete!(set, c_removed)
    end
  end

  return nothing
end

"""
    isdisjoint(c_a, c_b, is_sorted = false)

Checks whether two sets `c_a`, `c_b` have no common elements. Assumes by default that the sets are sorted.
"""
function isdisjoint(c_a::AbstractVector, c_b::AbstractVector; is_sorted = false)
  if !is_sorted
    sort!(c_a)
    sort!(c_b)
  end
  m = length(c_a)
  n = length(c_b)

  i = 1
  j = 1
  while i <= m && j <= n
    if c_a[i] < c_b[j]
      i += 1
    elseif c_b[j] < c_a[i]
      j += 1
    else
      return false
    end
  end
  return true
end


"""
    edge_metric(c_a::AbstractVector, c_b::AbstractVector, edge_weight::AbstractComplexityWeight)

Given two cliques `c_a` and `c_b` return a value for their edge weight.
"""
function edge_metric(c_a::T, c_b::T, edge_weight::AbstractComplexityWeight) where {T <: Union{AbstractVector, AbstractSet}}
  n_1 = length(c_a)
  n_2 = length(c_b)

  # merged block size
  n_m = union_dim(c_a, c_b)
  return compute_complexity_savings(n_1, n_2, n_m, edge_weight)
end

# Assuming the complexity of the projection is roughly O(n^3), how many operations are saved by projection the merged cliques
# instead of the individual cliques
compute_complexity_savings(n_1::Int, n_2::Int, n_m::Int, edge_weight::ComplexityWeight) = n_1^3 + n_2^3 - n_m^3


"""
    max_elem(A::SparseMatrixCSC)

Find the matrix indices (i, j) of the first maximum element among the elements stored in A.nzval
"""
function max_elem(A::SparseMatrixCSC)
  length(A.nzval) == 0 && throw(DomainError("Sparse matrix A doesn't contain any entries"))
  n = size(A, 2)

  ~, ind = findmax(A.nzval)
  row = A.rowval[ind]

  col = 0
  for c = 1:n
    col_indices = A.colptr[c]:A.colptr[c+1]-1
    if in(ind, col_indices)
      col = c
      break;
    end
  end
  return (row, col)
end


function edge_from_index(A::SparseMatrixCSC{Float64, Int}, ind::Int)
  # find the edge for that value
  row = A.rowval[ind]
  n = size(A, 2)
  col = 0
  for c = 1:n
    col_indices = A.colptr[c]:A.colptr[c+1]-1
    if in(ind, col_indices)
      col = c
      break;
    end
  end
  return (row, col)
end



" Given two cliques `c1` and `c2` in the tree `t`, return the parent clique first."
function determine_parent(t::SuperNodeTree, c1::Int, c2::Int)
  if in(c2, t.snd_child[c1])
    return c1, c2
  else
    return c2, c1
  end
end

"""
find_neighbors(edges::SparseMatrixCSC, c::Int)

Find all the cliques connected to `c` which are given by the nonzeros in `(c, 1:c-1)` and `(c+1:n, c)`.
"""
function find_neighbors(edges::SparseMatrixCSC, c::Int)
  neighbors = zeros(Int, 0)
  m, n = size(edges)
  # find all nonzero columns in row c up to column c
  if c > 1
   neighbors = vcat(neighbors, findall(x -> x!= 0, edges[c, 1:c-1]))
 end
 # find all nonzero entries in column c below c
 if c < n
  rows = edges.rowval[edges.colptr[c]:edges.colptr[c+1]-1]
  if edges.colptr[c] <= edges.colptr[c+1] - 1
    neighbors = vcat(neighbors, rows)
  end
end
return neighbors
end

function clique_intersections!(E::SparseMatrixCSC, snd::Array{Set{Int}, 1})
  # iterate over the nonzeros of the connectivity matrix E which represents the clique graph and replace the value by
  # |C_i ∩ C_j|
  rows = rowvals(E)
  for col in 1:size(E, 2)
    for j in nzrange(E, col)
      row = rows[j]
      E[row, col] = intersect_dim(snd[row], snd[col])
    end
  end
  return nothing
end


"""
    kruskal!(E::SparseMatrixCSC, num_cliques::Int)

Kruskal's algorithm to find a maximum weight spanning tree from the clique intersection graph.

 `E[i,j]` holds the cardinalities of the intersection between two cliques (i, j). Changes the entries in the connectivity matrix `E` to a negative
 value if an edge between two cliques is included in the max spanning tree.

 This is a modified version of https://github.com/JuliaGraphs/LightGraphs.jl/blob/master/src/spanningtrees/kruskal.jl
 """
function kruskal!(E::SparseMatrixCSC, num_cliques::Int)
  num_initial_cliques = size(E, 2)
  connected_c = DataStructures.IntDisjointSets(num_initial_cliques)

  I, J, V = findnz(E)
  # sort the weights and edges from maximum to minimum value
  p = sortperm(V, rev = true)
  I = I[p]
  J = J[p]
  V = V[p]
  num_edges_found = 0
  # iterate through edges (I -- J) with decreasing weight
  for k = 1:length(V)
    row = I[k]
    col = J[k]
    if !in_same_set(connected_c, row, col)
      union!(connected_c, row, col)
      # we indicate an edge in the MST with a positive value in E (all other values are >= 0)
      E[row, col] = -1.0
      num_edges_found += 1
      # break when all cliques are connected in one tree
      num_edges_found >= num_cliques - 1 && break
    end
  end
  return nothing
end


function assign_children!(snd_par::Array{Int, 1}, snd_child::Array{Set{Int}, 1}, c::Int, edges::SparseMatrixCSC)
  # determine neighbors
  neighbors = find_neighbors(edges, c)
  for n in neighbors
    # conditions that there is a edge in the MST and that n is not the parent of c
    if edges[max(c, n), min(c, n)] == -1.0 && snd_par[c] != n
      snd_par[n] = c
      push!(snd_child[c], n)
      assign_children!(snd_par, snd_child, n, edges)
    end
  end
  return nothing
end


" Given the maximum weight spanning tree represented by `E`, determine a parent structure `snd_par` for the clique tree."
function determine_parent_cliques!(snd_par::Array{Int, 1}, snd_child::Array{Set{Int}, 1}, cliques::Array{Set{Int}, 1}, post::Array{Int, 1}, E::SparseMatrixCSC)
  # vertex with highest order
  v = post[end]
  c = 0
  # Find clique that contains that vertex
  for (k, clique) in enumerate(cliques)
    if v ∈ clique
      # set that clique to the root
      snd_par[k] = 0
      c = k
      break
    end
  end

  # recursively assign children to cliques along the MST defined by E
  assign_children!(snd_par, snd_child, c, E)
  return nothing
end

" Traverse the clique tree in descending topological order and split the clique sets into supernodes and separators."
function split_cliques!(snd::Array{Set{Int}, 1}, sep::Array{Set{Int}, 1}, snd_par::Array{Int,1}, snd_post::Array{Int}, num_cliques::Int)

  # travese in topological decending order through the clique tree and split the clique in snd and sep
  for j = 1:1:(num_cliques - 1)
    c_ind = snd_post[j]
    par_ind = snd_par[c_ind]

    # find intersection of clique with parent
    # FIXME: Can be made a lot faster by using non-standard functions
    sep[c_ind] = intersect(snd[c_ind], snd[par_ind])
    snd[c_ind] = filter!(x -> x ∉ sep[c_ind], snd[c_ind])
  end
  return nothing
end


"""
    clique_tree_from_graph!(clique_graph::SuperNodeTree)

Given the cliques and edges of a clique graph, this function computes a valid clique tree.

This is necessary to perform the psd completion step of the dual variable after solving the problem.
"""
function clique_tree_from_graph!(t::SuperNodeTree)
  # a clique tree is a maximum weight spanning tree of the clique graph where the edge weight is the cardinality of the intersection between two cliques
  # compute intersection value for each edge in clique graph
  # clique_intersections!(t.strategy.edges, t.snd)
  clique_intersections!(t.strategy.edges, t.snd)

  # find a maximum weight spanning tree of the clique graph using Kruskal's algorithm
  kruskal!(t.strategy.edges, t.num)

  # determine the root clique of the clique tree (it can be any clique, but we AbstractComplexityWeightose the clique that contains the vertex with the highest order)
  determine_parent_cliques!(t.snd_par, t.snd_child, t.snd, t.post, t.strategy.edges)

  # recompute a postorder for the supernodes
  t.snd_post = post_order(t.snd_par, t.snd_child, t.num)

  t.sep = [Set{Int}() for i = 1:length(t.snd)]
  # split clique sets back into seperators and supernodes
  split_cliques!(t.snd, t.sep, t.snd_par, t.snd_post, t.num)
  t.strategy.clique_tree_recomputed = true
  return nothing
end


"Return the number of elements in s ∩ s2."
function intersect_dim(s::Set, s2::Set)
  if length(s) < length(s2)
        sa = s
        sb = s2
    else
        sa = s2
        sb = s
    end
    dim = 0
    for e in sa
        e in sb && (dim += 1)
    end
    return dim
end



" Find the size of the set `A ∪ B` under the assumption that `A` and `B` only have unique elements."
function union_dim(A::T, B::T)  where {T <: Union{AbstractVector, AbstractSet}}
  dim = length(A)
  for elem in B
    !in(elem, A) && (dim += 1)
  end
  return dim
end


"Compute the amount of fill-in created by merging two cliques with the respective supernode and separator dimensions."
function fill_in(dim_clique_snd::Int, dim_clique_sep::Int, dim_par_snd::Int, dim_par_sep::Int)
  dim_par = dim_par_snd + dim_par_sep
  dim_clique = dim_clique_snd + dim_clique_sep
  return ((dim_par - dim_clique_sep) * (dim_clique - dim_clique_sep))::Int
end

max_snd_size(dim_clique_snd::Int, dim_par_snd::Int) = max(dim_clique_snd, dim_par_snd)
clique_dim(t, c_ind) = length(t.snd[c_ind]), length(t.sep[c_ind])


"Store information about the merge of the two merge candidates `cand`."
function log_merge!(t::SuperNodeTree, do_merge::Bool, cand::Array{Int, 1})
  t.merge_log.clique_pairs = vcat(t.merge_log.clique_pairs, cand')
  push!(t.merge_log.decisions, do_merge)
  do_merge && (t.merge_log.num += 1)
  return nothing
end


"""
    print_merge_logs(ws; print_max = 10)

Print the merge logs for each decomposed cone in the problem.
"""
function print_merge_logs(ws; print_max = 10)
  sp_arr = ws.ci.sp_arr
  # print the merge logs for each sparsity pattern
  println(">>> Merge Logs:")
  for (iii, sp) in enumerate(sp_arr)
    m_log = sp.sntree.merge_log
    println("Sparsity Pattern Nr. $(iii), Graph Size: $(length(sp.sntree.par))")
    println("\t Num merges: $(m_log.num)\n\t Num decisions: $(length(m_log.decisions))\n\t Pairs:")

    # only print a certain number
    print_max = min(print_max, length(m_log.decisions))
    for jjj = 1:print_max
      p = m_log.clique_pairs[jjj, :]
      println("\t \tC_A: $(p[1]), C_B: $(p[2]) ($(m_log.decisions[jjj]))")
    end

  end
end

"Return an array of `MergeLog`s corresponding to the decomposed cones."
get_merge_logs(ws) = map(x -> x.sntree.merge_log, ws.ci.sp_arr)
