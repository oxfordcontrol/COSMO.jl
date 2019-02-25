
mutable struct SuperNodeTree
	snd::Array{Int64,1} #vertices of supernodes stored in one array
	snptr::Array{Int64,1} # vertices of supernode i are stored in snd[snptr[i]:snptr[i+1]-1]
	snd_par::Array{Int64,1}  # parent of supernode k is supernode j=snd_par[k]
	snd_post::Array{Int64,1} # post order of supernodal elimination tree
	post::Array{Int64} # post ordering of the vertices in elim tree σ(j) = v
	par::Array{Int64}
	cliques::Array{Int64,1} #vertices of cliques stored in one array
	chptr::Array{Int64,1} #points to the indizes where new clique starts in cliques array
	nBlk::Array{Int64,1} #sizes of submatrizes defined by each clique
	function SuperNodeTree(L)
		par = etree(L)
		child = child_from_par(par)
		post = post_order(par, child)
		# sorting of children from highest order one to lowest make sure that a node always gets
		# added to the supernode with the highest rep vertex
		sort_children!(child, post)
		# faster algorithm according to Vandenberghe p.308
		degrees = higher_degrees(L)
		snd, snptr, snd_par = find_supernodes(par, child, post, degrees)

		# TODO: amalgamate supernodes
		snd_child = child_from_par(snd_par)
		# # a post-ordering of the elimination tree is needed to make sure that in the
		# # psd completion step the matrix is completed from lower-right to top-left
		 snd_post = post_order(snd_par, snd_child)

		cliques, chptr, nBlk = find_cliques(L, snd, snptr, snd_par, post)

		new(snd, snptr, snd_par, snd_post, post, par, cliques, chptr, nBlk)
	end

end

function sort_children!(child, post)
	for v = 1:length(child)
		child[v] = sort(child[v], by = x -> post[x], rev = true)
	end
end
# -------------------------------------
# FUNCTION DEFINITIONS
# -------------------------------------
# given v=σ^-1(i) it returns i=σ(v)
function invert_order(sigma::Array{Int64,1})
	sigma_inv = zeros(Int64, length(sigma))
	for iii=1:length(sigma)
		sigma_inv[sigma[iii]] = iii
	end
	return sigma_inv
end

function num_cliques(sntree::SuperNodeTree)
	return length(sntree.chptr)
end

# elimination tree algorithm from H.Liu - A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# function etree_liu(g)
# 	N = length(g.adjacency_list)
# 	par = zeros(Int64,N)
# 	ancestor = zeros(Int64,N)

# 	elemSequence = g.reverse_ordering[collect(1:N)]
# 	for iii in elemSequence
# 		for vk in g.adjacency_list[iii]
# 			if g.ordering[vk] < g.ordering[iii]
# 				r = vk
# 				while (ancestor[r] != 0) && (ancestor[r] != iii)
# 					t = ancestor[r]
# 					ancestor[r] = iii
# 					r = t
# 				end
# 				if ancestor[r] == 0
# 					ancestor[r] = iii
# 					par[r] = iii
# 				end
# 			end
# 		end
# 	end

# 	return par
# end

# simplified version of my own elimination tree algorithm with simplified data structure (fastest)
function etree(L)
	N = size(L, 1)
	par = zeros(Int64, N)
	# loop over Vertices of graph
	for i=1:N
		value = i
		# number of i-neighbors with order higher than order of node i
		par_ = find_parent_direct(L, i)
		par[i] = par_
	end
	return par
end

# perform a depth-first-search to determine the post order of the tree defined by parent and children vectors
function post_order(par::Array{Int64,1}, child::Array{Array{Int64,1}})
	N = length(par)
	order = zeros(Int64, N)
	root = findall(x ->x == 0, par)[1]
	stack = Array{Int64}(undef, 0)
	iii = N
	push!(stack, root)
	while !isempty(stack)
		v = pop!(stack)
		order[v] = iii
		iii -= 1
		push!(stack, child[v]...)
	end
	post = collect(1:N)
	sort!(post, by = x-> order[x])
	return post
end



function child_from_par(par::Array{Int64,1})

	child = [Array{Int64}(undef, 0) for i = 1:length(par)]
	for i = 1:length(par)
		par_ = par[i]
		par_ != 0 && push!(child[par_], i)
	end
	return child
end

function get_snd(sntree::SuperNodeTree, ind::Int64)
	N = length(sntree.snptr)
	if ind == N
		return sntree.snd[sntree.snptr[ind]:end]
	else
		return sntree.snd[sntree.snptr[ind]:sntree.snptr[ind + 1] - 1]
	end
end

function get_clique(sntree::SuperNodeTree, ind::Int64)
	N = length(sntree.chptr)
	if ind == N
		return sntree.cliques[sntree.chptr[ind]:end]
	else
		return sntree.cliques[sntree.chptr[ind]:sntree.chptr[ind + 1] - 1]
	end
end

function print_cliques(sntree::SuperNodeTree)
	N = length(sntree.chptr)
	chptr = sntree.chptr
	println("Cliques of Graph:")
	for iii = 1:N
		if iii != N
			clique = sntree.cliques[chptr[iii]:chptr[iii + 1] - 1]
		else
			clique = sntree.cliques[chptr[iii]:end]
		end
		println("$(iii): $(clique)")
	end
end

function print_supernodes(sntree::SuperNodeTree)
	N = length(sntree.snptr)
	snptr = sntree.snptr
	println("Supernodes of Graph:")
	for iii = 1:N
		if iii != N
			sn = sntree.snd[snptr[iii]:snptr[iii + 1] - 1]
		else
			sn = sntree.snd[snptr[iii]:end]
		end
		println("$(iii): $(sn)")
	end
end


function check_degree_condition(v::Int64, w::Int64, degrees::Array{Int64,1})
	return degrees[v] > degrees[w] - 1
end


# Algorithm from A. Poten and C. Sun: Compact Clique Tree Data Structures in Sparse Matrix Factorizations (1989)
function pothen_sun(par::Array{Int64,1}, child::Array{Array{Int64,1}}, post::Array{Int64,1}, degrees::Array{Int64,1})
	N = length(par)
	sn_ind = -1 * ones(Int64, N) # if snInd[v] < 0 then v is a rep vertex, otherwise v ∈ supernode[snInd[v]]
	supernode_par = -1 * ones(Int64, N)

	# go through vertices in post_order
	for v in post
		child_ind = 0
		# check degree condition for all of v's childs from highest to lowest

		for (iii, w) in enumerate(child[v])
			# if not deg+(v) > deg+(w) - 1 for a certain w, set u to be w in snd(u), add v to snd(u)
			if !check_degree_condition(v, w, degrees)
				sn_ind[w] < 0 ? (u = w) : (u = sn_ind[w])
				sn_ind[v] = u
				child_ind = iii
				break
			end
		end

		# if v is still a rep vertex (i.e. above loop didnt find a child that fulfilled degree condition)
		if sn_ind[v] < 0
			u = v
		end

		for (iii, w) in enumerate(child[v])
			if iii != child_ind
				sn_ind[w] < 0 ? (x = w) : x = sn_ind[w]
				supernode_par[x] = u
			end
		end
	end

	# representative vertices
	reprv = findall(x-> x < 0, sn_ind)
	# vertices that are the parent of representative vertices
	repr_par = supernode_par[reprv]
	# take into account that all non-representative arrays are removed from the parent structure
	sn_par = zeros(Int64, length(reprv))

	for (iii, rp) in enumerate(repr_par)
		ind = findfirst(x -> x == rp, reprv)
		ind == nothing && (ind = 0)
		sn_par[iii] = ind
	end

	return sn_par, sn_ind
end

function find_supernodes(par::Array{Int64,1}, child::Array{Array{Int64,1}}, post::Array{Int64,1}, degrees::Array{Int64,1})
	supernode_par, snInd = pothen_sun(par, child, post, degrees)
	# number of vertices
	N = length(par)
	# number of representative vertices == number of supernodes
	Nrep = length(supernode_par)

	snode = zeros(Int64, N)
	snode_list = [Array{Int64}(undef, 0) for i = 1:N]
	snptr = zeros(Int64,Nrep)

	for iii in post
		f = snInd[iii]
		if f < 0
			push!(snode_list[iii], iii)

		else
			push!(snode_list[f], iii)
		end
	end

	kkk = 1
	jjj = 1
	for (iii, list) in enumerate(snode_list)
		len = length(list)
		if len > 0
			snptr[jjj] = kkk
			snode[kkk:kkk + len - 1] = list
			kkk += len
			jjj += 1
		end
	end
	return snode, snptr, supernode_par

end

function find_cliques(L, snodes::Array{Int64,1}, snptr::Array{Int64,1}, supernode_par::Array{Int64,1}, post::Array{Int64,1})
	postInv = invert_order(post)

	Nc = length(supernode_par)
	cliques = Array{Int64}(undef, 0)
	nBlk = zeros(Int64, Nc)
	chptr = zeros(Int64, Nc)
	jjj = 1

	for iii = 1:Nc
		if iii < Nc
			vRep = snodes[snptr[iii]:snptr[iii + 1] - 1][1]
		else
			vRep = snodes[snptr[iii]:end][1]
		end
		adjPlus = find_higher_order_neighbors(L, vRep)
		deg = length(adjPlus) + 1
		cliques = [cliques; vRep; adjPlus]
		chptr[iii] = jjj
		nBlk[iii] = Base.power_by_squaring(length(adjPlus) + 1, 2)
		jjj += deg
	end

	return cliques, chptr, nBlk

end

" Given a sparse lower triangular matrix L find the neighboring vertices of v"
function find_neighbors(L::SparseMatrixCSC, v::Int64)
	L = Symmetric(L)
	col_ptr = L.colptr
	row_val = L.rowval
	return row_val[col_ptr[v]:col_ptr[v + 1] - 1]
end

function find_higher_order_neighbors(L::SparseMatrixCSC, v::Int64)
	v == size(L, 1) && return 0
	col_ptr = L.colptr
	row_val = L.rowval
	return row_val[col_ptr[v]:col_ptr[v + 1] - 1]
end


function find_parent_direct(L::SparseMatrixCSC, v::Int64)
	v == size(L, 1) && return 0
	col_ptr = L.colptr
	row_val = L.rowval
	return row_val[col_ptr[v]]
end



# findall the cardinality of adj+(v) for all v in V
function higher_degrees(L::SparseMatrixCSC)
	N = size(L, 1)
	degrees = zeros(Int64, N)
	col_ptr = L.colptr
	row_val = L.rowval

	for v = 1:N-1
		degrees[v] = col_ptr[v + 1] - col_ptr[v]
	end
	return degrees
end


# -------------------------------------
# FUNCTION DEFINITIONS
# -------------------------------------
function find_graph!(ci, rows::Array{Int64, 1}, N::Int64, C::AbstractConvexSet)
	row_val, col_val = row_ind_to_matrix_indices(rows, N, C)
	F = QDLDL.qdldl(sparse(row_val, col_val, ones(length(row_val))), logical = true)
	ci.L = F.L
	return F.p
end


# given an array [rows] that represent the nonzero entries of a vectorized NxN matrix,
# return the rows and columns of the nonzero entries of the original matrix
function row_ind_to_matrix_indices(rows::Array{Int64,1}, N::Int64, ::PsdCone{Float64})
	row_val = zeros(Int64, length(rows))
	col_val = zeros(Int64, length(rows))
	for (ind, r) in enumerate(rows)
		_rem = mod(r, N)
		fl = fld(r, N)
		if _rem == 0
			row_val[ind] = N
			col_val[ind] = fl
		else
			row_val[ind] = _rem
			col_val[ind] = fl + 1
		end
	end
	return row_val, col_val
end

# given an array [rows] that represent the nonzero entries of the vectorized upper triangular part of a NxN matrix,
# return the rows and columns of the nonzero entries of the original matrix
function row_ind_to_matrix_indices(rows::Array{Int64,1}, N::Int64, ::COSMO.PsdConeTriangle{Float64})
	#  allocate conservative here since we don't know how many diagonal entries are contained in row
	row_val = zeros(Int64, 2 * length(rows) + N)
	col_val = zeros(Int64, 2 * length(rows) + N)
	ind = 1
	_step = 1
	for (iii, r) in enumerate(rows)
		i, j = COSMO.svec_to_mat(r)

		row_val[ind] = i
		col_val[ind] = j
		ind += 1

		if i != j
			row_val[ind] = j
			col_val[ind] = i
			ind += 1
		end
	end
	return row_val[1:ind - 1], col_val[1:ind - 1]
end

# Given a linear index find the col of the corresponding upper triangular matrix
function svec_get_col(x::Int64)
	c = (sqrt(8 * x + 1) - 1) / 2
	if c % 1 != 0
		return Int((floor(c) + 1))
	else
		return Int(c)
	end
end

# Given a linear index find the row of the corresponding upper triangular matrix
function svec_get_row(x::Int64)
	c = get_col(x) - 1
	k = (c + 1) * c / 2
	return (x - k)::Int64
end

# Given a linear index find the row and column of the corresponding upper triangular matrix
function svec_to_mat(ind::Int64)
	c = svec_get_col(ind)
	k = (c - 1) * c / 2
	r = Int(ind - k)
	return r::Int64, c::Int64
end






