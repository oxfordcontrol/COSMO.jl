function _contains(convex_sets::COSMO.CompositeConvexSet, t::Type{<:COSMO.AbstractConvexCone})
  for set in convex_sets.sets
    if typeof(set) <: t
      return true
    end
  end
  return false
end

function chordal_decomposition!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}
  # do nothing if no psd cones present in the problem
  if !_contains(ws.p.C, DecomposableCones)
    ws.ci.decompose = false
    return nothing
  end
  ws.ci = ChordalInfo{T}(ws.p, ws.settings)

  find_sparsity_patterns!(ws)

  if ws.ci.num_decomposable > 0

    # Do transformation similar to clique tree based transformation in SparseCoLo
    if ws.settings.compact_transformation
      augment_clique_based!(ws)
    else
      # find transformation matrix H and new composite convex set
      find_decomposition_matrix!(ws)
      # augment the original system
      augment_system!(ws)
    end
    pre_allocate_variables!(ws)
    ws.ci.decompose = true
    ws.states.IS_CHORDAL_DECOMPOSED = true
  else
    ws.ci.decompose = false
  end
  nothing
end

# analyse PSD cone constraints for chordal sparsity pattern
function find_sparsity_patterns!(ws::COSMO.Workspace)
  row_ranges = get_set_indices(ws.p.C.sets)
  sp_ind = 1
  for (k, C) in enumerate(ws.p.C.sets)
    psd_row_range = row_ranges[k]
    csp = find_aggregate_sparsity(ws.p.A, ws.p.b, psd_row_range, C)
    sp_ind = analyse_sparsity_pattern!(ws.ci, csp, ws.p.C.sets, C, k, psd_row_range, sp_ind, ws.settings.merge_strategy)
  end
end

analyse_sparsity_pattern!(ci::ChordalInfo, csp::Array{Int, 1}, sets::Vector{AbstractConvexSet}, C::AbstractConvexSet, k::Int, psd_row_range::UnitRange{Int}, sp_ind::Int, merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}}) = sp_ind

function analyse_sparsity_pattern!(ci::ChordalInfo, csp::Array{Int, 1}, sets::Vector{AbstractConvexSet}, C::DecomposableCones{T}, k::Int, psd_row_range::UnitRange{Int}, sp_ind::Int, merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}}) where {T <: Real}
  if length(csp) < C.dim
    return _analyse_sparsity_pattern(ci, csp, sets, C, k, psd_row_range, sp_ind, merge_strategy)
  else
   sets[k] = COSMO.DenseEquivalent(C, C.dim)
   return sp_ind
 end
end

function _analyse_sparsity_pattern(ci::ChordalInfo{T}, csp::Array{Int, 1}, sets::Vector{AbstractConvexSet}, C::Union{PsdCone{T}, PsdConeTriangle{T}}, k::Int, psd_row_range::UnitRange{Int}, sp_ind::Int, merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}}) where {T <: AbstractFloat}
  ordering, nz_ind_map = find_graph!(ci, csp, C.sqrt_dim, C)
  sp = COSMO.SparsityPattern(ci.L, C.sqrt_dim, ordering, merge_strategy, psd_row_range, k, nz_ind_map)
  # if after analysis of SparsityPattern & clique merging only one clique remains, don't bother decomposing
  if num_cliques(sp.sntree) == 1
    sets[k] = DenseEquivalent(C, C.dim)
    return sp_ind
  else
    ci.sp_arr[sp_ind] = sp
    push!(ci.psd_cones_ind, k)
    ci.num_decomposable += 1
    return sp_ind + 1
  end
end

DenseEquivalent(C::COSMO.PsdCone{T}, dim::Int) where {T <: AbstractFloat} = COSMO.DensePsdCone{T}(dim)
DenseEquivalent(C::COSMO.PsdConeTriangle{T, T}, dim::Int) where {T <: AbstractFloat} = COSMO.DensePsdConeTriangle{T, T}(dim)
DenseEquivalent(C::COSMO.PsdConeTriangle{T, Complex{T}}, dim::Int) where {T <: AbstractFloat} = COSMO.DensePsdConeTriangle{T, Complex{T}}(dim)

function nz_rows(a::SparseMatrixCSC{T}, ind::UnitRange{Int}, DROP_ZEROS_FLAG::Bool) where {T <: AbstractFloat}
  DROP_ZEROS_FLAG && dropzeros!(a)
  active = falses(length(ind))
  for r in a.rowval
    if in(r, ind)
      active[r - ind.start + 1] = true
    end
  end
  active
end

function number_of_overlaps_in_rows(A::SparseMatrixCSC{T}) where {T <: AbstractFloat}
  # sum the entries row-wise
  numOverlaps = sum(A, dims = 2)
  ri = findall(x -> x > 1, numOverlaps)
  return ri, numOverlaps[ri]
end


function find_aggregate_sparsity(A::SparseMatrixCSC{T}, b::AbstractVector{T}, ind::UnitRange{Int}, C::DecomposableCones{T}) where {T <: AbstractFloat}

  AInd_logical = nz_rows(A, ind, false)

  # explicitly flag all the terms corresonding to the cone diagonal
  for i = 1:C.sqrt_dim
    AInd_logical[vec_dim(i, C)] = true
  end

  AInd = findall(AInd_logical)

  # commonZeros = AInd[find(x->x==0,b[AInd])]
  bInd = findall(x -> x != 0, view(b, ind))
  commonNZeros = union(AInd, bInd)
  return commonNZeros
end
find_aggregate_sparsity(A::SparseMatrixCSC{T}, b::AbstractVector{T}, ind::UnitRange{Int}, C::AbstractConvexSet{T}) where {T <: AbstractFloat} = Int[]


"""
    reverse_decomposition!(ws::COSMO.Workspace, settings::COSMO.Settings)

After the problem has beend solved, undo the chordal decomposition.

Depending on the kind of transformation that was used, this involves:

- Reassembling the original matrix S from its blocks
- Reassembling the dual variable MU and performing a positive semidefinite completion.
"""
function reverse_decomposition!(ws::COSMO.Workspace{T}, settings::COSMO.Settings{T}) where {T <: AbstractFloat}

  mO = ws.ci.originalM
  nO = ws.ci.originalN
  vars = Variables{T}(mO, nO, ws.ci.originalC)
  vars.x .= ws.vars.x[1:nO]

  if settings.compact_transformation
    # reassemble the original variables s and μ
    add_sub_blocks!(vars.s, ws.vars.s, vars.μ, ws.vars.μ, ws.ci, ws.p.C, ws.ci.originalC, ws.ci.cone_map)
  else
    H = ws.ci.H
    vars.s  .= SplitVector{T}(H * ws.vars.s[mO + 1:end], ws.ci.originalC)
    fill_dual_variables!(ws, vars)
  end

  ws.p.C = ws.ci.originalC
  # if user requests, perform positive semidefinite completion on entries of μ that were not in the decomposed blocks
  ws.vars = vars
  settings.complete_dual && psd_completion!(ws)

  return nothing
end

function fill_dual_variables!(ws::COSMO.Workspace{T}, vars::COSMO.Variables{T}) where {T <: AbstractFloat}
  mO = ws.ci.originalM
  H = ws.ci.H

  # this performs the operation μ = sum H_k^T *  μ_k causing an addition of (identical valued) overlapping blocks
  vars.μ .= H * ws.vars.μ[mO + 1:end]

  # # to remove the overlaps we take the average of the values for each overlap by dividing by the number of blocks that overlap in a particular entry, i.e. number of 1s in each row of H
  rowInd, nnzs = number_of_overlaps_in_rows(H)

  for iii = 1:length(rowInd)
    ri = rowInd[iii]
    vars.μ[ri] = vars.μ[ri] / nnzs[iii]
  end
  return nothing
end

function add_sub_blocks!(s::SplitVector{T}, s_decomp::SplitVector{T}, μ::AbstractVector{T}, μ_decomp::AbstractVector{T}, ci::ChordalInfo{T}, C::CompositeConvexSet{T}, C0::CompositeConvexSet{T}, cone_map::Dict{Int, Int}) where {T <: AbstractFloat}
  sp_arr = ci.sp_arr
  row_start = 1 # the row pointer in the decomposed problem
  row_ranges = get_set_indices(C0.sets) # the row ranges of the same cone (or "parent" cone) in the original problem

  # iterate over all the cones of the decomposed problems and add the entries into the correct positions of the original problem
  for (k, C) in enumerate(C.sets)
    row_range = row_ranges[cone_map[k]]
    row_start = add_blocks!(s, μ, row_start, row_range, sp_arr, s_decomp, μ_decomp, C)
  end
  return nothing
end

function add_blocks!(s::SplitVector{T}, μ::AbstractVector{T}, row_start::Int, row_range::UnitRange{Int}, sp_arr::Array{SparsityPattern, 1}, s_decomp::SplitVector{T}, μ_decomp::AbstractVector{T}, C::AbstractConvexSet{T}) where {T <: AbstractFloat}

  @. s.data[row_range] = s_decomp.data[row_start:row_start + C.dim - 1]
  @. μ[row_range] = μ_decomp[row_start:row_start + C.dim - 1]
  return row_start + C.dim
end

function add_blocks!(s::SplitVector{T}, μ::AbstractVector{T}, row_start::Int, row_range::UnitRange{Int}, sp_arr::Array{SparsityPattern, 1}, s_decomp::SplitVector{T}, μ_decomp::AbstractVector{T}, C::DecomposableCones{T}) where {T <: AbstractFloat}
  # load the appropriate sparsity_pattern
  sp = sp_arr[C.tree_ind]
  sntree = sp.sntree
  ordering = sp.ordering
  #row_start = sp.row_range.start
  N = length(ordering)
  clique = map(v -> ordering[v], get_clique(sntree, C.clique_ind))
  sort!(clique)
  counter = 0
  for j in clique, i in clique
    if isa(C, PsdCone) || i <= j
      offset = COSMO.vectorized_ind(i, j, N, C) - 1
      s.data[row_range.start + offset] += s_decomp.data[row_start + counter]

      # notice: this overwrites the overlapping entries
      μ[row_range.start + offset] = μ_decomp[row_start + counter]
      counter += 1
    end
  end
  row_start += get_blk_rows(length(clique), C)

  return row_start
end

# The psd entries of μ that correspond to the zeros in s are not constrained by the problem
# however, in order to make the dual psd cone positive semidefinite we have to do a
# positive semidefinite completion routine to choose the values
function psd_completion!(ws::COSMO.Workspace)

  # loop over psd cones
  row_ranges = get_set_indices(ws.p.C.sets)
  sp_ind = 1
  for (kkk, C) in enumerate(ws.p.C.sets)
    sp_ind = complete!(ws.vars.μ, C, ws.ci.sp_arr, sp_ind, row_ranges[kkk])
  end

  return nothing
end

complete!(μ::AbstractVector{T}, ::AbstractConvexSet{T}, sp_arr::Array{SparsityPattern}, sp_ind::Int, rows::UnitRange{Int}) where {T <: AbstractFloat} = sp_ind

function complete!(μ::AbstractVector{T}, C::PsdCone{T}, sp_arr::Array{SparsityPattern}, sp_ind::Int, rows::UnitRange{Int}) where {T <: AbstractFloat}
  sp = sp_arr[sp_ind]

  μ_view = view(μ, rows)
  # make this y = -μ
  @. μ_view *= -one(T)

  M = reshape(μ_view, C.sqrt_dim, C.sqrt_dim)

  psd_complete!(M, C.sqrt_dim, sp.sntree, sp.ordering)

  @. μ_view *= -one(T)
  return sp_ind + 1
end

function complete!(μ::AbstractVector{T}, C::PsdConeTriangle{T}, sp_arr::Array{SparsityPattern}, sp_ind::Int, rows::UnitRange{Int}) where {T <: AbstractFloat}
  sp = sp_arr[sp_ind]

  μ_view = view(μ, rows)

  # I want to psd complete y, which is -μ
  populate_upper_triangle!(C.X, -μ_view, one(T) / sqrt(T(2)))
  psd_complete!(C.X, C.sqrt_dim, sp.sntree, sp.ordering)
  extract_upper_triangle!(C.X, μ_view, sqrt(2))
  @. μ_view *= -one(T)
  return sp_ind + 1
end


# positive semidefinite completion (from Vandenberghe - Chordal Graphs..., p. 362)
# input: A - positive definite completable matrix
function psd_complete!(A::AbstractMatrix{T}, N::Int, sntree::SuperNodeTree, p::Array{Int}) where {T <: AbstractFloat}

  # if a clique graph based merge strategy was used for this sparsity pattern, recompute a valid clique tree
  #recompute_clique_tree(sntree.strategy) && clique_tree_from_graph!(sntree, sntree.strategy)

  ip = invperm(p)

  As = Symmetric(A, :U)
  # permutate matrix based on ordering p (p must be a vector type), W is in the order that the cliques are based on
  W = copy(As[p, p])
  W = Matrix(W)
  num_cliques = COSMO.num_cliques(sntree)

  # go through supernode tree in descending order (given a post-ordering). This is ensured in the get_snd, get_sep functions
  for j = (num_cliques - 1):-1:1

    # in order to obtain ν, α the vertex numbers of the supernode are mapped to the new position of the permuted matrix
    # index set of snd(i) sorted using the numerical ordering i,i+1,...i+ni
    ν = get_snd(sntree, j)
    #clique = get_clique(sntree, snd_id)
    # index set containing the elements of col(i) \ snd(i) sorted using numerical ordering σ(i)
    α = get_sep(sntree, j)

    # index set containing the row indizes of the lower-triangular zeros in column i (i: representative index) sorted by σ(i)
    i = ν[1]
    η = collect(i+1:1:N)
    # filter out elements in lower triangular part of column i that are non-zero
    filter!(x -> !in(x, α) && !in(x, ν), η)

    Waa = W[α, α]
    Wαν = view(W, α, ν)
    Wηα = view(W, η, α)

    Y = zeros(length(α), length(ν))
    try
     Y[:, :] = Waa \ Wαν
   catch
    Waa_pinv = pinv(Waa)
    Y[:, :] = Waa_pinv * Wαν
  end

  W[η, ν] =  Wηα * Y
  # symmetry condition
  W[ν, η] = view(W, η, ν)'
end

# invert the permutation
A[:, :] =  W[ip, ip]
end
