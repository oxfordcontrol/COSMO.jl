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
`compute_complexity_savings(n_1::Int64, n_2::Int64, n_m::Int64, edge_weight::Subtype) --> Float64`

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

The (default) merge strategy based on the clique (intersection) graph ``\\mathcal{G}(\\mathcal{B}, \\xi)``, for a set of cliques ``\\mathcal{B} = \\{ \\mathcal{C}_1, \\dots, \\mathcal{C}_p\\}`` where the edge set ``\\xi`` is defined as
``\\xi = \\{ (\\mathcal{C}_{i},\\mathcal{C}_{j}) \\mid i \\neq j, \\; \\mathcal{C}_i, \\mathcal{C}_j \\in \\mathcal{B},\\; |\\mathcal{C}_i \\cap \\mathcal{C}_j| > 0 \\}``. In other words, we add an edge between every two cliques
 that overlap in at least one entry.

Moreover, given an edge weighting function ``e(\\mathcal{C}_i,\\mathcal{C}_j) = w_{ij}``, we compute a weight for each edge that quantifies the computational savings of merging the two cliques.
After the initial weights are computed, we merge cliques in a loop:

**while** clique graph contains positive weights:
- select two cliques with the highest weight ``w_{ij}``
- merge cliques ``\\rightarrow`` update clique graph
- recompute weights for updated clique graph

Custom edge weighting functions can be used by defining your own `CustomEdgeWeight <: AbstractEdgeWeight` and a corresponding `edge_metric` method. By default, the `ComplexityWeight <: AbstractEdgeWeight` is used which computes the weight based
on the cardinalities of the cliques: ``e(\\mathcal{C}_i,\\mathcal{C}_j)  = |\\mathcal{C}_i|^3 + |\\mathcal{C}_j|^3 - |\\mathcal{C}_i \\cup \\mathcal{C}_j|^3``.

See also: *Garstka, Cannon, Goulart - A clique graph based merging strategy for decomposable SDPs (2019)*
"""
mutable struct CliqueGraphMerge <: AbstractGraphBasedMerge
  stop::Bool
  edges::AbstractMatrix
  edge_weight::AbstractEdgeWeight
  clique_tree_recomputed::Bool
  function CliqueGraphMerge(; edge_weight = ComplexityWeight())
    new(false, spzeros(Float64, 0, 0), edge_weight, false)
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
  clique_ind::Int64
  t_fill::Int64
  t_size::Int64

  function ParentChildMerge(; t_fill = 8, t_size = 8)
    new(false, 2, t_fill, t_size)
  end
end

# Main clique merging routine:
# 1. initialise!() - tree and strategy
# 2. traverse() - find the next merge candidates
# 3. evalute() - evaluate wether to merge them or not
# 4. log() - log the decision
# 5. update!() - update the tree/graph and strategy
function _merge_cliques!(t::SuperNodeTree, strategy::AbstractMergeStrategy)
  initialise!(t, strategy)

  while !strategy.stop
    # find merge candidates
    cand = traverse(t, strategy)
    ordered_cand = copy(cand)
    # evaluate wether to merge the candidates
    do_merge = evaluate(t, strategy, cand)
    if do_merge
      # return cand in such a way that the remaining clique comes first
      ordered_cand[:, :] = merge_two_cliques!(t, cand, strategy)

    end
    log_merge!(t, do_merge, ordered_cand)

    # update strategy information after the merge
    update!(strategy, t, cand, ordered_cand, do_merge)
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
  t.snd_par = -ones(Int64, length(t.snd))

  # recomute a clique tree from the clique graph
  t.num > 1 && clique_tree_from_graph!(t)
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
function merge_child!(t::SuperNodeTree, cand::Array{Int64, 1})
  # determine which clique is the parent
  p, ch = determine_parent(t, cand[1], cand[2])

  # merge child's vertex sets into parent's vertex set
  push!(t.snd[p], t.snd[ch]...)
  t.snd[ch] = [] #necessary or just leave it
  t.sep[ch] = []

  # update parent structure
  @. t.snd_par[t.snd_child[ch]] = p
  t.snd_par[ch] = -1 #-1 instead of NaN, effectively remove that entry from the parent list

  # update children structure
  filter!(x -> x != ch, t.snd_child[p])
  push!(t.snd_child[p], t.snd_child[ch]...)
  t.snd_child[ch] = []

  # decrement number of cliques in tree
  t.num -= 1
  return [p; ch]
end
#
# # Merge two cliques that are in a sibling relationship
# function merge_sibling!(t::SuperNodeTree, cand::Array{Int64, 1})
#   c1 = cand[1]
#   c2 = cand[2]
#   # merge vertex set of cand[2] into cand[1]
#   push!(t.snd[c1], t.snd[c2]...)
#   t.snd[c2] = []
#   union!(t.sep[c1], t.sep[c2])
#   t.sep[c2] = []
#
#   @. t.snd_par[t.snd_child[c2]] = c1
#   t.snd_par[c2] = -1
#
#   # update children structure
#   push!(t.snd_child[c1], t.snd_child[c2]...)
#   t.snd_child[c2] = []
#
#   # decrement number of cliques in tree
#   t.num -= 1
#   return [c1; c2]
# end

"Given the clique graph `t` merge the two cliques with indices in `cand`."
function merge_two_cliques!(t::SuperNodeTree, cand::Array{Int64, 1}, strategy::AbstractGraphBasedMerge)
  c_1 = cand[1]
  c_2 = cand[2]
  # merge clique c_2 into c_1
  union!(t.snd[c_1], t.snd[c_2])
  sort!(t.snd[c_1])
  t.snd[c_2] = []

  # decrement number of cliques in tree
  t.num -= 1

  return [c_1; c_2]
end

"Given the clique tree `t` merge the two cliques with indices in `cand` as parent-child."
merge_two_cliques!(t::SuperNodeTree, cand::Array{Int64, 1}, strategy::ParentChildMerge) = merge_child!(t, cand)


# compute the edge set of the initial clique graph, only consider parent-child and sibling-sibling relationships
function initialise!(t::SuperNodeTree, strategy::CliqueGraphMerge)
  n_cliques = length(t.snd_par)
  strategy.edges = spzeros(Float64, n_cliques, n_cliques)

  # loop over all cliques
  for c = 1:n_cliques
    # brute force method for finding all edges to cliques that have some overlap
    compute_edges!(t, strategy, c)
  end
end

function initialise!(t::SuperNodeTree, strategy::ParentChildMerge)
  # start with node that has second highest order
  strategy.clique_ind = length(t.snd) - 1
end

"Find the next two cliques in the clique graph `t` to merge."
function traverse(t::SuperNodeTree, strategy::CliqueGraphMerge)
  # find maximum edge value in sparse edge matrix
  return max_elem(strategy.edges)
end

" Traverse tree `t` in descending topological order and return parent and clique (root has highest order)."
function traverse(t::SuperNodeTree, strategy::ParentChildMerge)
  c = t.snd_post[strategy.clique_ind]
  return [t.snd_par[c]; c]
end

"Decide whether to merge the two cliques with clique indices `cand`."
function evaluate(t, strategy::ParentChildMerge, cand::Array{Int64, 1})
  strategy.stop && return false
  par = cand[1]
  c = cand[2]
  dim_par_snd, dim_par_sep = COSMO.clique_dim(t, par)
  dim_clique_snd, dim_clique_sep = COSMO.clique_dim(t, c)
  return fill_in(dim_clique_snd, dim_clique_sep, dim_par_snd, dim_par_sep) <= strategy.t_fill || max_snd_size(dim_clique_snd, dim_par_snd) <= strategy.t_size
end

evaluate(t::SuperNodeTree, strategy::CliqueGraphMerge, cand::Array{Int64, 1}) = evaluate(t, strategy, cand, strategy.edge_weight)

function evaluate(t::SuperNodeTree, strategy::CliqueGraphMerge, cand::Array{Int64, 1}, edge_weight::AbstractComplexityWeight)
  do_merge = (strategy.edges[cand[1], cand[2]] >= 0)

  if !do_merge
    strategy.stop = true
  end
  return do_merge
end

" After a merge attempt, update the strategy information."
function update!(strategy::ParentChildMerge, t::SuperNodeTree, cand::Array{Int64, 1}, ordered_cand::Array{Int64, 1}, do_merge::Bool)
  # try to merge last node of order 1, then stop
  if strategy.clique_ind == 1
    strategy.stop = true
    # otherwise decrement node index
  else
    strategy.clique_ind -= 1
  end
end


function update!(strategy::CliqueGraphMerge, t::SuperNodeTree, cand::Array{Int64, 1}, ordered_cand::Array{Int64, 1}, do_merge::Bool)

  # After a merge operation update the information of the strategy
  if do_merge
    c_1_ind = cand[1]
    c_removed = cand[2]
    edges = strategy.edges
    n = size(edges, 2)

    c_1 = t.snd[c_1_ind]

    # recalculate edge values of all of c_1's and c_removed's neighbors and store in lower triangle
    n_arr = [find_neighbors(edges, c_1_ind); find_neighbors(edges, c_removed)]
    n_arr = unique!(n_arr)
    for n_ind in n_arr
      if n_ind != c_1_ind
        neigbhor = t.snd[n_ind]
        edges[max(c_1_ind, n_ind), min(c_1_ind, n_ind)] = COSMO.edge_metric(c_1, neigbhor, strategy.edge_weight)
      end
    end

    # remove entries in rows and columns that refer to edges that link to the removed clique
    edges[c_removed+1:n, c_removed] .= 0
    edges[c_removed, 1:c_removed] .= 0
    dropzeros!(edges)
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
    compute_edges!(t, strategy, c_a)

Computes the edge metric between clique `c_a` and all cliques that have some overlap with `c_a` and stores the result in `strategy.edges`.
"""
function compute_edges!(t::SuperNodeTree, strategy::CliqueGraphMerge, c_a_ind::Int64)
  # loop over all cliques (including empty ones and except c_a), and compute edge metric
  Nc = length(t.snd)
  c_a = t.snd[c_a_ind]
  for c_b_ind = c_a_ind+1:Nc
    c_b = t.snd[c_b_ind]
    if !isdisjoint(c_a, c_b; is_sorted = true)
      strategy.edges[c_b_ind, c_a_ind] = edge_metric(c_a, c_b, strategy.edge_weight)
    end
  end
end

"""
    edge_metric(c_a::AbstractVector, c_b::AbstractVector, edge_weight::AbstractComplexityWeight)

Given two cliques `c_a` and `c_b` return a value for their edge weight.
"""
function edge_metric(c_a::AbstractVector, c_b::AbstractVector, edge_weight::AbstractComplexityWeight)
  n_1 = length(c_a)
  n_2 = length(c_b)

  # merged block size
  n_m = union_dim(c_a, c_b)
  return compute_complexity_savings(n_1, n_2, n_m, edge_weight)
end

# Assuming the complexity of the projection is roughly O(n^3), how many operations are saved by projection the merged cliques
# instead of the individual cliques
compute_complexity_savings(n_1::Int64, n_2::Int64, n_m::Int64, edge_weight::ComplexityWeight) = n_1^3 + n_2^3 - n_m^3


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
  return [row; col]
end

" Given two cliques `c1` and `c2` in the tree `t`, return the parent clique first."
function determine_parent(t::SuperNodeTree, c1::Int64, c2::Int64)
  if in(c2, t.snd_child[c1])
    return c1, c2
  else
    return c2, c1
  end
end

"""
find_neighbors(edges::SparseMatrixCSC, c::Int64)

Find all the cliques connected to `c` which are given by the nonzeros in `(c, 1:c-1)` and `(c+1:n, c)`.
"""
function find_neighbors(edges::SparseMatrixCSC, c::Int64)
  neighbors = zeros(Int64, 0)
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

function clique_intersections!(E::SparseMatrixCSC, snd::Array{Array{Int64, 1}})
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
    kruskal!(E::SparseMatrixCSC, num_cliques::Int64)

Kruskal's algorithm to find a maximum weight spanning tree from the clique intersection graph.

 `E[i,j]` holds the cardinalities of the intersection between two cliques (i, j). Changes the entries in the connectivity matrix `E` to a negative
 value if an edge between two cliques is included in the max spanning tree.

 This is a modified version of https://github.com/JuliaGraphs/LightGraphs.jl/blob/master/src/spanningtrees/kruskal.jl
 """
function kruskal!(E::SparseMatrixCSC, num_cliques::Int64)
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

"Given a MWT represented by `edges`, find the children of clique `c` and store them in `snd_child."
function assign_children!(snd_par::Array{Int64, 1}, snd_child::Array{Array{Int64, 1}}, c::Int64, edges::SparseMatrixCSC)
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
function determine_parent_cliques!(snd_par::Array{Int64, 1}, snd_child::Array{Array{Int64, 1}}, cliques::Array{Array{Int64, 1}}, post::Array{Int64, 1}, E::SparseMatrixCSC)
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
function split_cliques!(snd::Array{Array{Int64,1}}, sep::Array{Array{Int64,1}}, snd_par::Array{Int64,1}, snd_post::Array{Int64}, num_cliques::Int64)

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
  clique_intersections!(t.strategy.edges, t.snd)

  # find a maximum weight spanning tree of the clique graph using Kruskal's algorithm
  kruskal!(t.strategy.edges, t.num)

  # determine the root clique of the clique tree (it can be any clique, but we AbstractComplexityWeightose the clique that contains the vertex with the highest order)
  determine_parent_cliques!(t.snd_par, t.snd_child, t.snd, t.post, t.strategy.edges)

  # recompute a postorder for the supernodes
  t.snd_post = post_order(t.snd_par, t.snd_child, t.num)

  # split clique sets back into seperators and supernodes
  split_cliques!(t.snd, t.sep, t.snd_par, t.snd_post, t.num)
  t.strategy.clique_tree_recomputed = true
  return nothing
end

" Find the size of the set `A ∩ B` under the assumption that `B` only has unique elements."
function intersect_dim(A::AbstractVector, B::AbstractVector)
  dim = 0
  for elem in B
    in(elem, A) && (dim += 1)
  end
  return dim
end

" Find the size of the set `A ∪ B` under the assumption that `A` and `B` only have unique elements."
function union_dim(A::AbstractVector, B::AbstractVector)
  dim = length(A)
  for elem in B
    !in(elem, A) && (dim += 1)
  end
  return dim
end


"Compute the amount of fill-in created by merging two cliques with the respective supernode and separator dimensions."
function fill_in(dim_clique_snd::Int64, dim_clique_sep::Int64, dim_par_snd::Int64, dim_par_sep::Int64)
  dim_par = dim_par_snd + dim_par_sep
  dim_clique = dim_clique_snd + dim_clique_sep
  return ((dim_par - dim_clique_sep) * (dim_clique - dim_clique_sep))::Int64
end

max_snd_size(dim_clique_snd::Int64, dim_par_snd::Int64) = max(dim_clique_snd, dim_par_snd)
clique_dim(t, c_ind) = length(t.snd[c_ind]), length(t.sep[c_ind])


"Store information about the merge of the two merge candidates `cand`."
function log_merge!(t::SuperNodeTree, do_merge::Bool, cand::Array{Int64, 1})
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
