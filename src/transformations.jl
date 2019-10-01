# -----------------------------------------
# Functions related to traditional decomposition
# -----------------------------------------
" Find the transformation matrix H to decompose the vector s into its parts and stacks them into sbar"
function find_decomposition_matrix!(ws)

  # allocate H and new decomposed cones
  n = COSMO.find_H_col_dimension(ws.p.C.sets, ws.ci.sp_arr)
  H_I = zeros(Int64, n)


  # find number of decomposed and total sets and allocate structure for new compositve convex set
  num_total, num_new_psd_cones = COSMO.num_cone_decomposition(ws)
  # decomposed_psd_cones = Array{COSMO.PsdCone}(undef, 0)
  C_new = Array{COSMO.AbstractConvexSet{Float64}}(undef, num_total)
  C_new[1] = COSMO.ZeroSet{Float64}(ws.ci.originalM)

  # loop over all convex sets and fill H and composite convex set accordingly
  row = 1
  entry = 1
  sp_ind = 1
  set_ind = 2
  for (kkk, C) in enumerate(ws.p.C.sets)
    set_ind, sp_ind, entry = COSMO.decompose!(H_I, C_new, set_ind, C, entry, row, ws.ci.sp_arr, sp_ind)
    row += C.dim
  end
  ws.p.C = COSMO.CompositeConvexSet(C_new)
  ws.ci.H = sparse(H_I, collect(1:n), ones(n))
end


function num_cone_decomposition(ws)
  num_sets = length(ws.p.C.sets)
  num_old_psd_cones = length(ws.ci.psd_cones_ind)
  num_new_psd_cones = 0
  for iii = 1:num_old_psd_cones
    sp = ws.ci.sp_arr[iii]
    num_new_psd_cones += COSMO.num_cliques(sp.sntree)
  end
  ws.ci.num_decom_psd_cones = ws.ci.num_psd_cones - ws.ci.num_decomposable + num_new_psd_cones
  # the total number is the number of original non-psd cones + number of decomposed psd cones + 1 zeroset for the augmented system
  num_total = num_sets - num_old_psd_cones + num_new_psd_cones + 1
  return num_total, num_new_psd_cones
end



function decompose!(H_I::Vector{Int64}, C_new, set_ind::Int64, C::COSMO.AbstractConvexSet, entry::Int64, row::Int64, sp_arr::Array{SparsityPattern}, sp_ind::Int64)
  #H[row:row + C.dim - 1, col:col + C.dim - 1] = sparse(1.0I, C.dim, C.dim)
  for i = 1:C.dim
    H_I[entry] = row + i - 1
    entry += 1
  end
  C_new[set_ind] = C

  return set_ind + 1, sp_ind, entry
end

# for a clique with nBlk entries, return the number of entries in the corresponding matrix
get_blk_rows(nBlk::Int64, C::COSMO.PsdCone{T}) where {T} = Base.power_by_squaring(nBlk, 2)
get_blk_rows(nBlk::Int64, C::COSMO.PsdConeTriangle{T}) where {T} = div((nBlk + 1) * nBlk, 2)

function decompose!(H_I::Vector{Int64}, C_new, set_ind::Int64,  C::DecomposableCones{ <: Real}, entry::Int64, row_start::Int64, sp_arr::Array{SparsityPattern}, sp_ind::Int64)
  sparsity_pattern = sp_arr[sp_ind]
  sntree = sparsity_pattern.sntree
  # original matrix size
  original_size = C.sqrt_dim

  for iii = 1:num_cliques(sntree)

    c = COSMO.get_clique(sntree, iii)
    # the graph and tree algorithms determined the cliques vertices of an AMD-permuted matrix. Since the location of the data hasn't changed in reality, we have to map the clique vertices back
    c = map(v -> sparsity_pattern.ordering[v], c)
    sort!(c)
    entry = COSMO.add_subblock_map!(H_I, c, original_size, entry, row_start, C)

    # create and add new cone for subblock
    num_rows = get_blk_rows(get_nBlk(sntree, iii), C)
    C_new[set_ind] = typeof(C)(num_rows)
    set_ind += 1
  end
  return set_ind, sp_ind + 1, entry
end



# fills the corresponding entries of H for clique c
function add_subblock_map!(H_I::Vector{Int64}, clique_vertices::Array{Int64}, original_size::Int64, entry::Int64, row_start::Int64, ::PsdCone{<: Real})
  for vj in clique_vertices
    for vi in clique_vertices
      row = mat_to_vec_ind(vi, vj, original_size)
      H_I[entry] = row_start + row - 1
      entry += 1
    end
  end
  return entry::Int64
end

function add_subblock_map!(H_I::Vector{Int64}, clique_vertices::Array{Int64}, original_size::Int64, entry::Int64,  row_start::Int64, ::PsdConeTriangle{<: Real})
  for vj in clique_vertices
    for vi in clique_vertices
      if vi <= vj
        row = mat_to_svec_ind(vi, vj)
        H_I[entry] = row_start + row - 1
        entry += 1
      end
    end
  end
  return entry::Int64
end


function find_H_col_dimension(sets, sp_arr)
  sum_cols = 0
  sp_arr_ind = 1
  for C in sets
    dim, sp_arr_ind = decomposed_dim(C, sp_arr, sp_arr_ind)
    sum_cols += dim
  end
  return sum_cols::Int64
end


function augment_system!(ws)
  _augment!(ws.p, ws.ci.H)
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



# -----------------------------------------
# Functions related to clique tree based transformation
# see: Kim: Exploiting sparsity in linear and nonlinear matrix inequalities via positive semidefinite matrix completion (2011), p.53
# -----------------------------------------

"""
    augment_clique_based!(ws::COSMO.Workspace{T})

If `settings.colo_transformation = true`, decompose the problem based on the clique tree.
"""
function augment_clique_based!(ws::COSMO.Workspace{T}) where {T}
  A = ws.p.A
  b = ws.p.b
  q = ws.p.q
  P = ws.p.P
  cones = ws.p.C.sets
  sp_arr = ws.ci.sp_arr


  # determine number of final augmented matrix and number of overlapping entries
  mA, nA, num_overlapping_entries = find_A_dimension(ws.p.model_size[2], ws.p.C.sets, ws.ci.sp_arr)
  # find number of decomposed and total sets and allocate structure for new compositve convex set
  num_total, num_new_psd_cones = COSMO.num_cone_decomposition(ws)

  C_new = Array{COSMO.AbstractConvexSet{Float64}}(undef, num_total - 1)

  # allocate memory for Aa_I, Aa_J, Aa_V, b_I
  nz = nnz(A)
  Aa_I = zeros(Int64, nz + 2 * num_overlapping_entries)
  Aa_J = extra_columns(nz + 2 * num_overlapping_entries, nz + 1, ws.p.model_size[2] + 1)
  Aa_V =  alternating_sequence(nz + 2 * num_overlapping_entries, nz + 1)
  findnz!(Aa_I, Aa_J, Aa_V, A)
  bs = sparse(b)
  b_I = zeros(Int64, length(bs.nzval))
  b_V = bs.nzval


  row_ranges = COSMO.get_set_indices(cones) # the row ranges of each cone in the original problem
  row_ptr = 1 # a continuously updated pointer to the start of the first entry of a cone in A_I
  overlap_ptr = nz + 1 # a continuously updated pointer to the next free entry to enter the rows for the +1, -1 overlap entries

  sp_ind = 1 # every decomposable cone is linked to a sparsity_pattern with index sp_ind
  set_ind = 1 # the set index counter set_ind is used to create a map of the set number in the original problem to the set number in the decomposed problem

  # Main Loop: Loop over the cones of the original problem and add entries to the row vectors Aa_I and b_I
  for (k, C) in enumerate(cones)
    row_range = row_ranges[k]
    row_ptr, overlap_ptr, set_ind, sp_ind = COSMO.add_entries!(Aa_I, b_I, C_new, row_ptr, A, bs, row_range, overlap_ptr, set_ind, sp_ind, sp_arr, C, k, ws.ci.cone_map)
  end

  ws.p.A = allocate_sparse_matrix(Aa_I, Aa_J, Aa_V, mA, nA)
  ws.p.b = Vector(SparseArrays._sparsevector!(b_I, b_V, mA))
  ws.p.P = blockdiag(P, spzeros(num_overlapping_entries, num_overlapping_entries))
  ws.p.q = vec([q; zeros(num_overlapping_entries)])
  ws.p.model_size[1] = size(ws.p.A, 1)
  ws.p.model_size[2] = size(ws.p.A, 2)
  ws.p.C = COSMO.CompositeConvexSet(C_new)
  return nothing
end

"Given the row, column, and nzval vectors and dimensions, assemble the sparse matrix `Aa` of the decomposed problem in a slightly more memory efficient way."
function allocate_sparse_matrix(Aa_I::Array{Int64, 1}, Aa_J::Array{Int64, 1}, Aa_V::Array{Float64, 1}, mA::Int64, nA::Int64)
  csrrowptr = zeros(Int64, mA + 1)
  csrcolval = zeros(Int64, length(Aa_I))
  csrnzval = zeros(Float64, length(Aa_I))
  klasttouch = zeros(Int64, nA + 1)
  csccolptr = zeros(Int64, nA + 1)
  # sort_col_wise!(Aa_I, Aa_V, A.colptr, size(A, 2))
  #Aa = SparseMatrixCSC{Float64, Int64}(mA, nA, Aa_J, Aa_I, Aa_V)
  Aa = SparseArrays.sparse!(Aa_I, Aa_J, Aa_V, mA, nA, +, klasttouch, csrrowptr, csrcolval, csrnzval, csccolptr, Aa_I, Aa_V )
end



"Return the dimension of the problem after a clique tree based decomposition, given the sparsity patterns in `sp_arr`."
function find_A_dimension(n_original::Int64, sets, sp_arr)
  num_cols = n_original
  num_overlapping_entries = 0
  num_rows = 0
  sp_arr_ind = 1
  for C in sets
    dim, overlaps, sp_arr_ind = decomposed_dim(C, sp_arr, sp_arr_ind)
    num_rows += dim
    num_overlapping_entries += overlaps
  end
  return num_rows::Int64, (num_cols + num_overlapping_entries)::Int64, num_overlapping_entries
end


# This method handles all non-decomposable sets C. For these non-decomposable sets you just have to shift the row indices of all entries by an offset
function add_entries!(A_I::Array{Int64, 1}, b_I::Array{Int64, 1}, C_new::Array{COSMO.AbstractConvexSet{Float64}, 1}, row_ptr::Int64, A0::SparseMatrixCSC, b0::SparseVector, row_range::UnitRange{Int64}, overlap_ptr::Int64, set_ind::Int64, sp_ind::Int64,
  sp_arr::Array{SparsityPattern, 1}, C::AbstractConvexSet, k::Int64, cone_map::Dict{Int64, Int64})

  m, n = size(A0)

  # let's handle vector b first
  offset = row_ptr - row_range.start

  row_range_col = COSMO.get_rows(b0, row_range)
  if row_range_col != nothing
    for k in row_range_col
      b_I[k] = b0.nzind[k] + offset
    end
  end

  # also create new set once
  C_new[set_ind] = C
  cone_map[set_ind] = k

  for col = 1:n
    ptr = copy(row_ptr)
    # indices which store the rows in column for C in A0
    row_range_col = COSMO.get_rows(A0, col, row_range)
    for k in row_range_col
      A_I[k] = A0.rowval[k] + offset
    end
  end


  return row_ptr + C.dim, overlap_ptr, set_ind + 1, sp_ind
end

# This method handles decomposable cones, e.g. PsdConeTriangle. The row vectors A_I and b_I have to be edited in such a way that entries of one clique are appear continously
function add_entries!(A_I::Array{Int64, 1}, b_I::Array{Int64, 1}, C_new::Array{COSMO.AbstractConvexSet{Float64}, 1}, row_ptr::Int64, A0::SparseMatrixCSC, b0::SparseVector, row_range::UnitRange{Int64}, overlap_ptr::Int64, set_ind::Int64, sp_ind::Int64,
  sp_arr::Array{SparsityPattern, 1},  C::PsdConeTriangle{Float64}, k::Int64,  cone_map::Dict{Int64, Int64})

  sp = sp_arr[sp_ind] # The SparsityPattern correspondig to this cone C
  sntree = sp.sntree # The Supernodal Elimination Tree that stores information about the cliques
  ordering = sp.ordering # A reordering that was applied by the LDL routine
  N_v = length(ordering)
  m, n = size(A0) # Dimensions of the original problem

  # determine the row ranges for each of the subblocks
  clique_to_rows = COSMO.clique_rows_map(row_ptr, sntree, C)

  # loop over cliques in descending topological order
  for iii = num_cliques(sntree):-1:1

    # get supernodes and seperators and undo the reordering
    sep = map(v -> ordering[v], get_sep(sntree, iii))
    isa(sep, Array{Any, 1}) && (sep = Int64[])
    snd = map(v -> ordering[v], get_snd(sntree, iii))
    # compute sorted block indices (i, j, flag) for this clique with an information flag whether an entry (i, j) is an overlap
    block_indices = COSMO.get_block_indices(snd, sep, N_v)

    # If we encounter an overlap with a parent clique we have to be able to find the location of the overlapping entry
    # Therefore load and reorder the parent clique
    if iii == num_cliques(sntree)
      par_clique = Int64[]
      par_rows = 0:0
    else
      par_ind = COSMO.get_clique_par(sntree, iii)
      par_rows = clique_to_rows[par_ind]
      par_clique = map(v -> ordering[v], get_clique_by_ind(sntree, par_ind))
      sort!(par_clique)
    end

    # Loop over all the columns and shift the rows in A_I and b_I according to the clique strucutre
    for col = 1:n
      row_range_col = COSMO.get_rows(A0, col, row_range)
      row_range_b = col == 1 ? COSMO.get_rows(b0, row_range) : 0:0
      overlap_ptr = add_clique_entries!(A_I, b_I, A0.rowval, b0.nzind, block_indices, par_clique, par_rows, col, C.sqrt_dim, row_ptr, overlap_ptr, row_range, row_range_col, row_range_b)
    end

    # create and add new cone for subblock
    num_rows = get_blk_rows(get_nBlk(sntree, iii), C)
    cone_map[set_ind] = k
    C_new[set_ind] = typeof(C)(num_rows, sp_ind, iii)
    row_ptr += num_rows
    set_ind += 1
  end
  return row_ptr, overlap_ptr, set_ind, sp_ind + 1
end

" Loop over all entries (i, j) in the clique and either set the correct row in `A_I` and `b_I` if (i, j) is not an overlap or add an overlap column with (-1 and +1) in the correct positions."
function add_clique_entries!(A_I::Array{Int64, 1}, b_I::Array{Int64, 1}, A_rowval::Array{Int64}, b_nzind::Array{Int64, 1}, block_indices::Array{Tuple{Int64, Int64, Int64},1},  par_clique::Array{Int64, 1}, par_rows::UnitRange{Int64}, col::Int64,  C_sqrt_dim::Int64, row_ptr::Int64, overlap_ptr::Int64, row_range::UnitRange{Int64}, row_range_col::UnitRange{Int64}, row_range_b::UnitRange{Int64})
  counter = 0
  for block_idx in block_indices
    new_row_val = row_ptr + counter
    # a block index that corresponds to an overlap
    if block_idx[3] == 0
      if col == 1
        i = block_idx[1]
        j = block_idx[2]
        A_I[overlap_ptr] = new_row_val # this creates the +1 entry
        A_I[overlap_ptr + 1] = par_rows.start + COSMO.parent_block_indices(par_clique, i, j) - 1 # this creates the -1 entry
        overlap_ptr += 2
      end
    else
      # (i, j) of the clique
      i = block_idx[1]
      j = block_idx[2]
      # k = svec(i, j)
      k = COSMO.mat_to_svec_ind(i, j)
      modify_clique_rows!(A_I, k, A_rowval, C_sqrt_dim, new_row_val, row_range, row_range_col)
      col == 1 && modify_clique_rows!(b_I, k, b_nzind, C_sqrt_dim, new_row_val, row_range, row_range_b)
    end
  counter += 1
  end
  return overlap_ptr
end

" Given the nominal entry position `k = svec(i, j)` find and modify with `new_row_val` the actual location of that entry in the global row vector `rowval`,"
function modify_clique_rows!(A_I::Array{Int64, 1}, k::Int64, rowval::Array{Int64, 1}, C_sqrt_dim::Int64, new_row_val::Int64, row_range::UnitRange{Int64}, row_range_col::UnitRange{Int64})
  row_0 = COSMO.get_row_index(k, rowval, C_sqrt_dim, row_range, row_range_col)
  # row_0 happens when (i, j) references an edge that was added by merging cliques, the corresponding value will be zero
  # and can be disregarded
  if row_0 != 0
    A_I[row_0] = new_row_val
  end
  return nothing
end

# function modify_clique_rows_b!(b_I::Array{Int64, 1}, k::Int64, b_nzind::Array{Int64, 1}, C_sqrt_dim::Int64, new_row_val::Int64, row_range::UnitRange{Int64}, row_range_b::UnitRange{Int64} )
#   row_0 = COSMO.get_row_index(k, b_nzind, C_sqrt_dim, row_range, row_range_b)
#   row_0 != 0 && (b_I[row_0] = new_row_val)
#   return nothing
# end


"Given the svec index `k` and an offset `row_range_col.start`, return the location of the (i, j)th entry in the row vector `rowval`."
function get_row_index(k::Int64, rowval::Array{Int64, 1}, sqrt_dim::Int64, row_range::UnitRange{Int64}, row_range_col::UnitRange{Int64})

  k_shift = row_range.start + k - 1

  # determine upper set boundary of where the row could be
  u = min(row_range_col.stop, row_range_col.start + k_shift - 1)

  # find index of first entry >= k, starting in the interval [l, u]
  # if no, entry is >= k, returns u + 1
  r = searchsortedfirst(rowval, k_shift, row_range_col.start, u, Base.Order.Forward)
  # if no r s.t. rowval[r] = k_shift was found that means that the (i, j)th entry represents an edded edge (zero) from clique merging
  if r > u || rowval[r] != k_shift
    return 0
  else
    return r
  end
end


# -----------------------------------------
# Utility Functions
# -----------------------------------------

" Find the index of k=svec(i, j) in the parent clique `par_clique`."
function parent_block_indices(par_clique::Array{Int64, 1}, i::Int64, j::Int64)
  ir = searchsortedfirst(par_clique, i)
  jr = searchsortedfirst(par_clique, j)
  return COSMO.mat_to_svec_ind(ir, jr)
end

"""
    get_block_indices(snd::Array{Int64}, sep::Array{Int64})

For a clique consisting of supernodes `snd` and seperators `sep`, compute all the indices (i, j) of the corresponding matrix block
in the format (i, j, flag) where flag is equal to 0 if entry (i, j) corresponds to an overlap of the clique and 1 otherwise.

`Nv` is the number of vertices in the graph that we are trying to decompose.
"""
function get_block_indices(snd::Array{Int64}, sep::Array{Int64}, Nv::Int64)
  N = length(sep) + length(snd)
  d = div(N * (N + 1), 2)

  block_indices = Array{Tuple{Int64, Int64, Int64}, 1}(undef, d)
  ind = 1

  for j in sep, i in sep
    if i <= j
      block_indices[ind] = (i, j, 0)
      ind += 1
    end
  end

  for j in snd, i in snd
    if i <= j
      block_indices[ind] = (i, j, 1)
      ind += 1
    end
  end

  for i in snd
    for j in sep
      block_indices[ind] = (min(i, j), max(i, j), 1)
      ind += 1
    end
  end

  sort!(block_indices, by = x -> x[2] * Nv + x[1] )
  return block_indices
end

decomposed_dim(C::AbstractConvexSet, sp_arr::Array{SparsityPattern}, sp_arr_ind::Int64) = (C.dim, 0, sp_arr_ind)
function decomposed_dim(C::DecomposableCones{ <: Real}, sp_arr::Array{SparsityPattern}, sp_arr_ind::Int64)
  sntree = sp_arr[sp_arr_ind].sntree
  dim, overlaps = get_decomposed_dim(sntree, C)
  return dim::Int64, overlaps::Int64, (sp_arr_ind + 1)::Int64
end

"Return the row ranges of each clique after the decomposition of `C` shifted by `row_start`."
function clique_rows_map(row_start::Int64, sntree::SuperNodeTree, C::DecomposableCones{<:Real})
  Nc = num_cliques(sntree)
  rows = Array{UnitRange{Int64}}(undef,  Nc)
  ind = zeros(Int64, Nc)
  for iii = Nc:-1:1
    num_rows = COSMO.get_blk_rows(COSMO.get_nBlk(sntree, iii), C)
    rows[iii] = row_start:row_start+num_rows-1
    ind[iii] = sntree.snd_post[iii]
    row_start += num_rows
  end
  return Dict(ind .=> rows)
end


function get_rows(b::SparseVector, row_range::UnitRange{Int64})
  rows = b.nzind
  if length(rows) > 0
    s = searchsortedfirst(rows, row_range.start)
    if rows[s] > row_range.stop || s == 0
        return nothing, 0:0
    else
      e = searchsortedlast(rows, row_range.stop)
      return s:e
    end
  else
    return nothing
  end

end

function get_rows(A::SparseMatrixCSC, col::Int64, row_range::UnitRange{Int64})
  colrange = A.colptr[col]:(A.colptr[col + 1]-1)

  # if the column has entries
  if colrange.start <= colrange.stop
    # create a view into the row values of column col
    rows = view(A.rowval, colrange)
    # find the rows within row_start:row_start+C.dim-1
    # s: index of first entry in rows >= row_start
    s = searchsortedfirst(rows, row_range.start)
    if rows[s] > row_range.stop || s == 0
      return 0:0
    else
      # e: index of last value in rows <= row_start + C.dim - 1
      e = searchsortedlast(rows, row_range.stop)
      return colrange[s]:colrange[e]
    end
  else
    return 0:0
  end
end

"Returns the appropriate amount of memory for `A.nzval`, including, starting from `n_start`, the (+1 -1) entries for the overlaps."
function alternating_sequence(total_length::Int64, n_start::Int64)
  v = ones(Float64, total_length)
  for i=n_start + 1:2:length(v)
    v[i] = -1
  end
  return v
end

"Returns the appropriate amount of memory for the columns of the augmented problem matrix `A`, including, starting from `n_start`, the columns for the (+1 -1) entries for the overlaps."
function extra_columns(total_length::Int64, n_start::Int64, start_val::Int64)
  v = zeros(Int64, total_length)
  for i = n_start:2:length(v)-1
    v[i] = start_val
    v[i + 1] = start_val
    start_val += 1
  end
  return v
end

function augmented_col_ptr(colptr::Array{Int64}, num_overlapping_entries::Int64)
  l = length(colptr)
  v = zeros(Int64, l + num_overlapping_entries)
  col = colptr[end] + 2
  for i = l+1:length(v)
    v[i] = col
    col += 2
  end
  return v
end

"Given a sparse matrix `S`, write the columns and non-zero values into the first `numnz` entries of `J` and `V`."
function findnz!(I::Vector{Ti}, J::Vector{Ti}, V::Vector{Tv}, S::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    numnz = nnz(S)
    count = 1
    @inbounds for col = 1 : S.n, k = S.colptr[col] : (S.colptr[col+1]-1)
       # I[count] = S.rowval[k]
        J[count] = col
        V[count] = S.nzval[k]
        count += 1
    end
end

function add_nz_entries!(J::Vector{Ti}, V::Vector{Tv}, S::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
  for i = 1:length(S.colptr)
    J[i] = S.colptr[i]
  end
  for i = 1:length(S.nzval)
    V[i] = S.nzval[i]
  end
  return nothing
end

"Given the row vector and colptr of a sparse matrix, sort the  entries within `rowvals` and `nzvals` withing the first `n0` columns."
function sort_col_wise!(rowval::Array{Int64, 1}, nzval::Array{Float64, 1}, colptr::Array{Int64, 1}, n0::Int64)
  for col = 1:n0
    r = colptr[col]:colptr[col + 1]-1
    row_view = view(rowval, r)
    nzval_view = view(nzval, r)
    p = sortperm(row_view)
    if !issorted(p)
      permute!(row_view, p)
      permute!(nzval_view, p)
    end
  end
  return nothing
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

# Converts the matrix element (i, j) of A ∈ m x n into the corresponding linear index of v = vec(A)
function mat_to_vec_ind(i::Int64, j::Int64, m::Int64)
  (i > m || i <= 0 || j <= 0) && throw(BoundsError("Indices outside matrix bounds."))
  return (j - 1) * m + i
end

# Converts the matrix element (i, j) of A ∈ m x n into the corresponding linear index of v = svec(A, ::UpperTriangle)
function mat_to_svec_ind(i::Int64, j::Int64)
  if i <= j
    return div((j - 1) * j, 2) + i
  else
    return div((i - 1) * i, 2) + j
  end
end

vectorized_ind(i::Int64, j::Int64, m::Int64, C::PsdCone{T}) where {T} = mat_to_vec_ind(i, j, m)
vectorized_ind(i::Int64, j::Int64, m::Int64, C::PsdConeTriangle{T}) where {T} = mat_to_svec_ind(i, j)

function svec_to_mat_ind(k::Int64)
  j = isqrt(2 * k)
  i = k - div((j - 1) * j, 2)
  return i, j
end

