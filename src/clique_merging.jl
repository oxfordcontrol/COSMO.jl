# -------------------------------------
# Functions related to clique merging
# -------------------------------------

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

get_merge_logs(ws) = map(x -> x.sntree.merge_log, ws.ci.sp_arr)

"""
NoMerge <: AbstractMergeStrategy

A strategy that does not merge cliques.
"""
struct NoMerge <: AbstractTreeBasedMerge end

abstract type AbstractEdgeScore end


struct RelIntersect <: AbstractEdgeScore end
abstract type AbstractComplexityScore <: AbstractEdgeScore end
struct ComplexityScore <: AbstractComplexityScore end


merge_cliques!(t::SuperNodeTree, strategy::NoMerge) = nothing


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

    t.num == 1 && break
    strategy.stop && break
    # update strategy information after the merge
    update!(strategy, t, cand, ordered_cand)
  end
return nothing
end

function merge_cliques!(t, strategy::AbstractTreeBasedMerge)
  # call main merge routine
  _merge_cliques!(t, strategy)
  # do some merge strategy specific post-processing of the tree
  Nc = t.num
  # the merging creates empty supernodes and seperators, recalculate a post order for the supernodes
  snd_post = post_order(t.snd_par, t.snd_child, Nc)
  t.snd_post = snd_post
  t.connectivity = spzeros(0, 0)
  return nothing
end

function merge_cliques!(t, strategy::AbstractGraphBasedMerge)
  # call main merge routine
  _merge_cliques!(t, strategy)
  # since for now we have a graph, not a tree a post ordering does not make sense. Therefore just number
  # the non-empty supernodes in t.snd
  t.snd_post = findall(x -> !isempty(x), t.snd)
  return nothing
end

# calculate block sizes (notice: the block sizes are stored in post order)
function calculate_block_dimensions!(t::SuperNodeTree, strategy::AbstractTreeBasedMerge)
  Nc = t.num
  t.nBlk = zeros(Nc)
  for iii = 1:Nc
    c = t.snd_post[iii]
    t.nBlk[iii] = Base.power_by_squaring(length(t.sep[c]) + length(t.snd[c]), 2)
  end
end


function calculate_block_dimensions!(t::SuperNodeTree, strategy::AbstractGraphBasedMerge)
  Nc = t.num
  t.nBlk = zeros(Nc)
  for iii = 1:Nc
    c = t.snd_post[iii]
    t.nBlk[iii] = Base.power_by_squaring(length(t.snd[c]), 2)
  end
end

function log_merge!(t::SuperNodeTree, do_merge::Bool, cand::Array{Int64, 1})
  t.merge_log.clique_pairs = vcat(t.merge_log.clique_pairs, cand')
  push!(t.merge_log.decisions, do_merge)
  do_merge && (t.merge_log.num += 1)
  return nothing
end

# Merge two cliques that are in a parent - child relationship
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

# Merge two cliques that are in a sibling relationship
function merge_sibling!(t::SuperNodeTree, cand::Array{Int64, 1})
  c1 = cand[1]
  c2 = cand[2]
  # merge vertex set of cand[2] into cand[1]
  push!(t.snd[c1], t.snd[c2]...)
  t.snd[c2] = []
  union!(t.sep[c1], t.sep[c2])
  t.sep[c2] = []

  @. t.snd_par[t.snd_child[c2]] = c1
  t.snd_par[c2] = -1

  # update children structure
  push!(t.snd_child[c1], t.snd_child[c2]...)
  t.snd_child[c2] = []

  # decrement number of cliques in tree
  t.num -= 1
  return [c1; c2]
end


"""
PairwiseMerge <: AbstractGraphBasedMerge

A merge strategy that calculates the edge metric `A ∩ B / A ∪ B` for every two cliques that are in a parent-child or sibling relationship. The resulting clique
graph is traversed from the highest edge metric to the lowest.
"""
mutable struct PairwiseMerge <: AbstractGraphBasedMerge
  stop::Bool
  edges::AbstractMatrix
  edge_score::AbstractEdgeScore
  recompute_ct_on_demand::Bool # if true only recompute clique tree from clique graph when psd_complete is called, i.e. dual variable has to be completed
  function PairwiseMerge(; edge_score = RelIntersect())
    new(false, spzeros(Float64, 0, 0), edge_score, true)
  end
end


"""
TreeTraversalMerge(t_fill = 8, t_size = 8) <: AbstractTreeBasedMerge

The merge strategy suggested in Sun / Andersen - Decomposition in conic optimization with partially separable structure (2013).
The initial clique tree is traversed in topological order and cliques are greedily merged to their parent if evaluate() returns true.
"""
mutable struct TreeTraversalMerge <: AbstractTreeBasedMerge
  stop::Bool
  clique_ind::Int64
  t_fill::Int64
  t_size::Int64

  function TreeTraversalMerge(; t_fill = 8, t_size = 8)
    new(false, 2, t_fill, t_size)
  end
end

"""
NakataMerge(zeta = 0.4) <: AbstractTreeBasedMerge

The merge strategy suggested in Nakata - Exploiting sparsity in semidefinite programming via matrix completion II (2001)
The initial clique tree is traversed in a depth first search and for each clique a conditition is checked if the clique's
 - any pair of children should be merged (sibling-merge)
 - the clique should be merged with a child (parent-child merge)
 - the clique and two children should be merged (3-family merge)
"""
mutable struct NakataMerge <: AbstractTreeBasedMerge
  stop::Bool
  clique_ind::Int64
  zeta::Int64


  function NakataMerge(; zeta = 0.4)
    new(false, 1, zeta)
  end
end





# compute the edge set of the initial clique graph, only consider parent-child and sibling-sibling relationships
function initialise!(t, strategy::PairwiseMerge)
  n_cliques = length(t.snd_par)
  strategy.edges = spzeros(Float64, n_cliques, n_cliques)

  # loop over all cliques
  for c = 1:n_cliques
    # brute force method for finding all edges to cliques that have some overlap
    compute_edges!(t, strategy, c)
  end
end

"""
isdisjoint(c_a, c_b)

Checks whether two sets c_a, c_b have no common elements.
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

Computes the edge metric between clique c_a and all cliques that have some overlap with c_a and stores the result in strategy.edges.
"""
function compute_edges!(t, strategy::PairwiseMerge, c_a_ind::Int64)
  # loop over all cliques (including empty ones and except c_a), and compute edge metric
  Nc = length(t.snd)
  c_a = t.snd[c_a_ind]
  for c_b_ind = c_a_ind+1:Nc
    c_b = t.snd[c_b_ind]
    if !isdisjoint(c_a, c_b; is_sorted = true)
      strategy.edges[c_b_ind, c_a_ind] = edge_metric(c_a, c_b, strategy.edge_score)
    end
  end
end

function edge_metric(c_a::AbstractVector, c_b::AbstractVector, edge_score::AbstractComplexityScore)
  n_1 = length(c_a)
  n_2 = length(c_b)

  # merged block size
  n_m = union_dim(c_a, c_b)
  return compute_complexity_savings(n_1, n_2, n_m, edge_score)
end

function edge_metric(c_a::AbstractVector, c_b::AbstractVector, edge_score::RelIntersect)
  return intersect_dim(c_a, c_b) / union_dim(c_a, c_b)
end


function traverse(t, strategy::PairwiseMerge)
  # find maximum edge value in sparse edge matrix
  return max_elem(strategy.edges)
end


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

merged_block_size(t::SuperNodeTree, c1::Int64, c2::Int64) = c1 < c2 ? block_size_child(t, c1, c2) : block_size_sibling(t, c1, c2)

# Given two cliques c1 and c2, return the parent clique first
function determine_parent(t::SuperNodeTree, c1::Int64, c2::Int64)
  if in(c2, t.snd_child[c1])
    return c1, c2
  else
    return c2, c1
  end
end

function block_size_child(t::SuperNodeTree, c1::Int64, c2::Int64)
  p, ch = determine_parent(t, c1, c2)
  return length(t.snd[p]) + length(t.sep[p]) + length(t.snd[ch])
end

function block_size_sibling(t::SuperNodeTree, c1::Int64, c2::Int64)
  return length(t.snd[c1]) + length(t.snd[c2]) + union_dim(t.sep[c1], t.sep[c2])
end

evaluate(t, strategy::PairwiseMerge, cand) = evaluate(t, strategy, cand, strategy.edge_score)


function evaluate(t, strategy::PairwiseMerge, cand, edge_score::RelIntersect)

  c1_ind = cand[1]
  c2_ind = cand[2]

  # case if only one clique remains
  if c1_ind == c2_ind
    strategy.stop = true
    return false
  end


  n_ops_diff = edge_metric(t.snd[c1_ind], t.snd[c2_ind], ComplexityScore())
  do_merge = (n_ops_diff >= 0)

  if !do_merge
    strategy.stop = true
  end
  return do_merge
end

function evaluate(t, strategy::PairwiseMerge, cand, edge_score::AbstractComplexityScore)
  do_merge = (strategy.edges[cand[1], cand[2]] >= 0)

  if !do_merge
    strategy.stop = true
  end
  return do_merge
end

# Assuming the complexity of the projection is roughly O(n^3), how many operations are saved by projection the merged cliques
# instead of the individual cliques
compute_complexity_savings(n_1::Int64, n_2::Int64, n_m::Int64, edge_score::ComplexityScore) = n_1^3 + n_2^3 - n_m^3

# Approximates the number of operations for one projection of all the cliques in the tree
compute_complexity(t::COSMO.SuperNodeTree) = sum(map(x -> x^3, t.nBlk))


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

function merge_two_cliques!(t::SuperNodeTree, cand::Array{Int64, 1}, strategy::AbstractTreeBasedMerge)
  # parent - child relationships are stored in upper triangle
  cand[1] < cand[2] ? merge_child!(t, cand) : merge_sibling!(t, cand)
end

# TreeTraversal merge strategy always merges parent-child pairs
merge_two_cliques!(t::SuperNodeTree, cand::Array{Int64, 1}, strategy::TreeTraversalMerge) = merge_child!(t, cand)

"""
find_neighbors(edges::SparseMatrixCSC, c::Int64)

Find all the cliques connected to c which are given by the nonzeros in (c, 1:c-1) and (c+1:n, c).
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


# After a merge operation update the information of the strategy
function update!(strategy::PairwiseMerge, t, cand, ordered_cand)
  c_1_ind = cand[1]
  c_removed = cand[2]
  edges = strategy.edges
  n = size(edges, 2)


  # ... and all to the removed clique
  edges[c_removed+1:n, c_removed] .= 0
  edges[c_removed, 1:c_removed] .= 0
  dropzeros!(strategy.edges)

  c_1 = t.snd[c_1_ind]
  # recalculate edge values of all of c_1's neighbors and store in lower triangle

  n_arr = find_neighbors(edges, c_1_ind)
  # println("c1: $(c_1_ind):")
  # @show(n_arr)
  for n_ind in n_arr
    neigbhor = t.snd[n_ind]
    edges[max(c_1_ind, n_ind), min(c_1_ind, n_ind)] = edge_metric(c_1, neigbhor, strategy.edge_score)
  end

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


# Kruskal's algorithm to find a maximum weight spanning tree from the clique intersection graph, E represents the cardinality of the
# intersection between two cliques
# Changes the entries in the connectivity matrix E to negative value if an edge between two cliques is included in the max spanning tree
# modified version of:https://github.com/JuliaGraphs/LightGraphs.jl/blob/master/src/spanningtrees/kruskal.jl
function kruskal!(E::SparseMatrixCSC)


  num_c = size(E, 2)
  connected_c = DataStructures.IntDisjointSets(num_c)

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
      num_edges_found >= num_c - 1 && break
    end
  end
  return nothing
end

  function assign_children!(snd_par::Array{Int64, 1}, snd_child::Array{Array{Int64, 1}}, c::Int64, edges::SparseMatrixCSC)
    # determine neighbors
    neighbors = find_neighbors(edges, c)
    for n in neighbors
      if n < c
        row = c
        col = n
      else
        row = n
        col = c
      end

      # conditions that there is a edge in the MST and that n is not the parent of c
      if edges[row, col] == -1.0 && snd_par[c] != n
        snd_par[n] = c
        push!(snd_child[c], n)
        assign_children!(snd_par, snd_child, n, edges)
      end
    end
    return nothing
end


function determine_parent_cliques!(snd_par::Array{Int64, 1}, snd_child::Array{Array{Int64, 1}}, cliques::Array{Array{Int64, 1}}, post::Array{Int64, 1}, E::SparseMatrixCSC)
  # vertex with highest order
  v = post[end]
  c = 0
  # find clique that contains that vertex
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

recompute_clique_tree(strategy::AbstractTreeBasedMerge) = false
recompute_clique_tree(strategy::AbstractGraphBasedMerge) = strategy.recompute_ct_on_demand

"""
    clique_tree_from_graph!(tree)

Given the cliques and edges of a clique graph, this function computes a valid clique tree. This is necessary to perform the psd completion step.
"""
function clique_tree_from_graph!(t::SuperNodeTree)
  # a clique tree is a maximum weight spanning tree of the clique graph where the edge weight is the cardinality of the intersection between two cliques
  # compute intersection value for each edge in clique graph
  clique_intersections!(t.strategy.edges, t.snd)

  # find a maximum weight spanning tree of the clique graph using Kruskal's algorithm
  kruskal!(t.strategy.edges)

  # determine the root clique of the clique tree (the clique that contains the vertex with the highest order)
  determine_parent_cliques!(t.snd_par, t.snd_child, t.snd, t.post, t.strategy.edges)

  # recompute a postorder
  t.snd_post = post_order(t.snd_par, t.snd_child)

  # split clique sets back into seperators and supernodes
  split_cliques!(t.snd, t.sep, t.snd_par, t.snd_post, t.num)
  return nothing
end

clique_tree_from_graph!(t::SuperNodeTree, strategy::AbstractTreeBasedMerge) = nothing

# # computes the edge metric for cliques C_a and C_b: (C_a ∩ C_b) / (C_a ∪ C_b) where C_a is the parent of C_b
# function edge_metric_parent(res, sep, c_a, c_b, edge_score::RelIntersect)
#   return length(sep[c_b]) / (length(sep[c_a]) + length(res[c_a]) + length(res[c_b]))
# end

# # computes the edge metric for cliques C_a and C_b: (C_a ∩ C_b) / (C_a ∪ C_b) where C_a is the sibling of C_b
# function edge_metric_siblings(res, sep, c_a, c_b, edge_score::RelIntersect)
#   return intersect_dim(sep[c_a], sep[c_b]) / (length(res[c_a]) + length(res[c_b]) + union_dim(sep[c_a], sep[c_b]))
# end

# # computes the edge metric for cliques C_a and C_b in terms of how the number of projection operationss change when merged
# function edge_metric_parent(res, sep, c_a, c_b, edge_score::ComplexityScore)

#   n_1 = length(res[c_a]) + length(sep[c_a])
#   n_2 = length(res[c_b]) + length(sep[c_b])
#   # merged block size
#   n_m = n_1 + length(res[c_b])
#   return compute_complexity_savings(n_1, n_2, n_m)
# end

# # computes the edge metric for cliques C_a and C_b: (C_a ∩ C_b) / (C_a ∪ C_b) where C_a is the sibling of C_b
# function edge_metric_siblings(res, sep, c_a, c_b, edge_score::ComplexityScore)
#   n_1 = length(res[c_a]) + length(sep[c_a])
#   n_2 = length(res[c_b]) + length(sep[c_b])
#   # merged block size
#   n_m = length(res[c_a]) + length(res[c_b]) + union_dim(sep[c_a], sep[c_b])
#   return compute_complexity_savings(n_1, n_2, n_m)
# end

# Finds the size of the set A ∩ B under the assumption that B only has unique elements
function intersect_dim(A::AbstractVector, B::AbstractVector)
  dim = 0
  for elem in B
    in(elem, A) && (dim += 1)
  end
  return dim
end

# Finds the size of the set A ∪ B under the assumption that A and B only have unique elements
function union_dim(A::AbstractVector, B::AbstractVector)
  dim = length(A)
  for elem in B
    !in(elem, A) && (dim += 1)
  end
  return dim
end

function initialise!(t, strategy::TreeTraversalMerge)
  # start with node that has second highest order
  strategy.clique_ind = length(t.snd) - 1
end

function initialise!(t, strategy::NakataMerge)
  # start with root (clique with highest order)
  strategy.clique_ind = length(t.snd)
end

# traverse tree in descending topological order and return parent and clique, root has highest order
function traverse(t, strategy::TreeTraversalMerge)
  c = t.snd_post[strategy.clique_ind]
  return [t.snd_par[c]; c]
end

function fill_in(dim_clique_snd::Int64, dim_clique_sep::Int64, dim_par_snd::Int64, dim_par_sep::Int64)
  dim_par = dim_par_snd + dim_par_sep
  dim_clique = dim_clique_snd + dim_clique_sep
  return ((dim_par - dim_clique_sep) * (dim_clique - dim_clique_sep))::Int64
end
max_snd_size(dim_clique_snd::Int64, dim_par_snd::Int64) = max(dim_clique_snd, dim_par_snd)

clique_dim(t, c_ind) = length(t.snd[c_ind]), length(t.sep[c_ind])

function evaluate(t, strategy::TreeTraversalMerge, cand)
  strategy.stop && return false
  par = cand[1]
  c = cand[2]
  dim_par_snd, dim_par_sep = COSMO.clique_dim(t, par)
  dim_clique_snd, dim_clique_sep = COSMO.clique_dim(t, c)
  #println("par: $(par), c: $(c), fill_in :$(fill_in(dim_clique_snd, dim_clique_sep, dim_par_snd, dim_par_sep)), size: $( max_snd_size(dim_clique_snd, dim_par_snd))")
  return fill_in(dim_clique_snd, dim_clique_sep, dim_par_snd, dim_par_sep) <= strategy.t_fill || max_snd_size(dim_clique_snd, dim_par_snd) <= strategy.t_size
end

function evaluate(t, strategy::NakataMerge, cand)
  strategy.stop && return false
  if length(cand) == 3
    return _evaluate(t, c1, c2, c3)
  else
    return _evaluate(t, c1, c2)
  end
end

function _evaluate(t, strategy, c_r, c_s)
  zeta = strategy.zeta
  intersecct_dim = intersect_dim(c_r, c_s)
  c_r_dim = t.snd[c_r] + t.sep[c_r]
  c_s_dim = t.snd[c_s] + t.sep[c_s]
  return min(intersect_dim / c_r_dim, intersect_dim / c_s_dim) >= zeta
end

function _evaluate(t, strategy, c_q, c_r, c_s)
  zeta = strategy.zeta
  intersecct_dim = intersect_dim(c_a, c_b)
  c_a_dim = t.snd[c_a] + t.sep[c_a]
  c_b_dim = t.snd[c_b] + t.sep[c_b]
  return min(intersect_dim / c_a_dim, intersect_dim / c_b_dim) >= zeta
end


function update!(strategy::TreeTraversalMerge, t, cand, ordered_cand)
  # try to merge last node of order 1, then stop
  if strategy.clique_ind == 1
    strategy.stop = true
    # otherwise decrement node index
  else
    strategy.clique_ind -= 1
  end
end