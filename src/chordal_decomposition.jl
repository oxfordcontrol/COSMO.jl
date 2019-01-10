function _contains(convex_sets::COSMO.CompositeConvexSet, t::Type{<:COSMO.AbstractConvexCone})
  for set in convex_sets.sets
    if typeof(set) == t
      return true
    end
  end
  return false
end

function chordal_decomposition!(ws::COSMO.Workspace)
  ws.ci = ChordalInfo{Float64}(ws.p)
  settings = ws.settings
  problem = ws.p

  # do nothing if no psd cones present in the problem
  if !_contains(problem.C, PsdCone{Float64})
    settings.decompose = false
    return nothing
  end

  find_sparsity_patterns!(ws)

  # find transformation matrix H and new composite convex set
  find_decomposition_matrix!(ws)

  # augment the original system
  augment_system!(ws)

  nothing
end

# analyses PSD cone constraints for chordal sparsity pattern
function find_sparsity_patterns!(ws)
  # find sparsity pattern, graphs, and clique sets for each cone
  for (iii, ind) in enumerate(ws.ci.psd_cones_ind)
    csp = find_aggregate_sparsity(ws.p.A, ws.p.b, ind)
    cDim = Int(sqrt(ind.stop - ind.start + 1))
    ws.ci.sp_arr[iii] = COSMO.SparsityPattern(csp, cDim, true)
  end
end

function nz_rows(a::SparseMatrixCSC, ind::UnitRange{Int64}, DROP_ZEROS_FLAG::Bool)
    DROP_ZEROS_FLAG && dropzeros!(a)
    active = falses(length(ind))
    for r in a.rowval
      if in(r, ind)
        active[r - ind.start + 1] = true
      end
    end
    return findall(active)
end

function number_of_overlaps_in_rows(A::SparseMatrixCSC)
  # sum the entries row-wise
  numOverlaps = sum(A, dims=2)
  ri = findall(x -> x > 1, numOverlaps)
  return ri, numOverlaps[ri]
end

function find_aggregate_sparsity(A, b, ind)
  AInd = nz_rows(A, ind, false)
  # commonZeros = AInd[find(x->x==0,b[AInd])]
  bInd = findall(x -> x != 0, view(b, ind))
  commonNZeros = union(AInd, bInd)
  return commonNZeros
end

function vec_to_mat_ind(ind::Int64, n::Int64)
  ind > n^2 && error("Index ind out of range.")
  ind == 1 && (return 1, 1)

  r = ind % n

  if r == 0
    j = Int(ind / n)
    i = n
  else
    j = Int(floor(ind / n) + 1)
    i = r
  end
  return i, j
end

# Converts the matrix element (i,j) of A ∈ m x n into the corresponding row value of v = vec(A)
function mat_to_vec_ind(i::Int64, j::Int64, m::Int64)
  (i > m || i <= 0 || j <= 0) && throw(BoundsError("Indices outside matrix bounds."))
  return (j - 1) * m + i
end


# function finds the transformation matrix H to decompose the vector s into its parts and stacks them into sbar
function find_decomposition_matrix!(ws)

  # allocate H and new decomposed cones
  m, n = COSMO.find_H_dimension(ws.ci)
  ws.ci.H = spzeros(m, n)

  # find number of decomposed and total sets and allocate structure for new compositve convex set
  num_total, num_new_psd_cones = COSMO.num_cone_decomposition(ws)
  # decomposed_psd_cones = Array{COSMO.PsdCone}(undef, 0)
  C_new = Array{COSMO.AbstractConvexSet{Float64}}(undef, num_total)
  C_new[1] = COSMO.ZeroSet{Float64}(m)

  # loop over all convex sets and fill H and composite convex set accordingly
  row = 1
  col = 1
  sp_ind = 1
  set_ind = 2
  for (kkk, C) in enumerate(ws.p.C.sets)
    set_ind, sp_ind, col = COSMO.decompose!(ws.ci.H, C_new, set_ind, C, row, col, ws.ci.sp_arr, sp_ind)
    row += C.dim
  end
  ws.p.C = COSMO.CompositeConvexSet(C_new)
end

function decompose!(H::SparseMatrixCSC, C_new, set_ind::Int64, C::COSMO.AbstractConvexSet, row::Int64, col::Int64, sp_arr::Array{SparsityPattern}, sp_ind::Int64)
  H[row:row + C.dim - 1, col:col + C.dim - 1] = sparse(1.0I, C.dim, C.dim)
  C_new[set_ind] = C

  return set_ind + 1, sp_ind, col + C.dim
end

function decompose!(H::SparseMatrixCSC, C_new, set_ind::Int64,  C::COSMO.PsdCone{Float64}, row_start::Int64, col_start::Int64, sp_arr::Array{SparsityPattern}, sp_ind::Int64)
  sntree = sp_arr[sp_ind].sntree
  # original matrix size
  original_size = C.sqrt_dim

  for iii = 1:num_cliques(sntree)
    # new stacked size
    block_size = Int(sqrt(sntree.nBlk[iii]))

    c = COSMO.get_clique(sntree, iii)
    sort!(c)
    col_start = COSMO.add_subblock_map!(H, c, block_size, original_size, row_start, col_start)

    # create and add new cone for subblock
    C_new[set_ind] = COSMO.PsdCone(block_size^2)
    set_ind += 1
  end
  return set_ind, sp_ind + 1, col_start
end

# fills the corresponding entries of H for clique c
function add_subblock_map!(H::SparseMatrixCSC, clique_vertices::Array{Int64}, block_size::Int64, original_size::Int64, row_shift::Int64, col_shift::Int64)
  num = 1
  for (iii, vj) in enumerate(clique_vertices)
    for (jjj, vi) in enumerate(clique_vertices)
      row = mat_to_vec_ind(vi, vj, original_size)
      H[row_shift + row - 1, col_shift + num - 1] = 1.
      num += 1
    end
  end
  return (col_shift + block_size^2)::Int64
end

function find_H_dimension(ci)
  m = ci.originalM
  num_decomp_cones = length(ci.psd_cones_ind)
  stacked_sizes = zeros(Int64, num_decomp_cones)
  for iii = 1:num_decomp_cones
    stacked_sizes[iii] = sum(ci.sp_arr[iii].sntree.nBlk)
  end

  num_other_rows = ci.originalM - sum(map(x -> length(x), ci.psd_cones_ind))

  # length of stacked vector sbar
  n = num_other_rows + sum(stacked_sizes)
  return m, n
end

function num_cone_decomposition(ws)
  num_sets = length(ws.p.C.sets)
  num_old_psd_cones = length(ws.ci.psd_cones_ind)
  num_new_psd_cones = 0
  for (iii, sp) in enumerate(ws.ci.sp_arr)
    num_new_psd_cones += COSMO.num_cliques(sp.sntree)
  end
  # the total number is the number of original non-psd cones + number of decomposed psd cones + 1 zeroset for the augmented system
  num_total = num_sets - num_old_psd_cones + num_new_psd_cones + 1
  return num_total, num_new_psd_cones
end

function augment_system!(ws)
  _augment!(ws.p, ws.ci.H)

  # increase the variable dimension
  ws.vars = Variables{Float64}(ws.p.model_size[1], ws.p.model_size[2], ws.p.C)
  nothing
end

function _augment!(problem, H::SparseMatrixCSC)
  mH, nH = size(H)
  m, n = size(problem.A)
  problem.P = blockdiag(problem.P, spzeros(nH, nH))
  problem.q = vec([problem.q; zeros(nH)])
  problem.A = [problem.A H; spzeros(nH, n) -sparse(1.0I, nH, nH)]
  problem.b = vec([problem.b; zeros(nH)])
  problem.model_size[1] = size(problem.A, 1)
  problem.model_size[2] = size(problem.A, 2)
  nothing
end


function reverse_decomposition!(ws::COSMO.Workspace, settings::COSMO.Settings)
  mO = ws.ci.originalM
  nO = ws.ci.originalN

  H = ws.ci.H

  vars = Variables{Float64}(mO, nO, ws.ci.originalC)
  vars.x .= ws.vars.x[1:nO]
  vars.s  .= SplitVector{Float64}(H * ws.vars.s[mO + 1:end], ws.ci.originalC)

  # fill dual variables such that μ_k  = H_k μ for k=1,...,p
  fill_dual_variables!(ws, vars)

  # if user requests, perform positive semidefinite completion on entries of μ that were not in the decomposed blocks
  settings.complete_dual && psd_completion!(ws, vars)

  ws.vars = vars
  ws.p.C = ws.ci.originalC
  return nothing
end

function fill_dual_variables!(ws::COSMO.Workspace,vars::Variables)
  mO = ws.ci.originalM
  H = ws.ci.H

  # this performs the operation μ = sum H_k^T *  μ_k causing an addition of (identical valued) overlapping blocks
  vars.μ .= H * ws.vars.μ[mO + 1:end]

  # to remove the overlaps we take the average of the values for each overlap by dividing by the number of blocks that overlap in a particular entry, i.e. number of 1s in each row of H
  rowInd,nnzs = number_of_overlaps_in_rows(H)

  for iii=1:length(rowInd)
    ri = rowInd[iii]
    vars.μ[ri] .= ws.vars.μ[ri] / nnzs[iii]
  end
end

# complete the dual variable
function psd_completion!(ws::COSMO.Workspace)
  return nothing
end
