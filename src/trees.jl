
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
	function SuperNodeTree(g::Graph)
		par = etree(g)
		child = child_from_par(par)
		post = post_order(par, child)
		# faster algorithm according to Vandenberghe p.308
		degrees = higher_degrees(g)

		snd, snptr, snd_par = find_supernodes(par, child, post, degrees)

		# TODO: amalgamate supernodes
		snd_child = child_from_par(snd_par)
		snd_post = post_order(snd_par, snd_child)

		cliques, chptr, nBlk = find_cliques(g, snd, snptr, snd_par, post)

		new(snd, snptr, snd_par, snd_post, post, par, cliques, chptr, nBlk)
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
function etree_liu(g::Graph)
	N = length(g.adjacency_list)
	par = zeros(Int64,N)
	ancestor = zeros(Int64,N)

	elemSequence = g.reverse_ordering[collect(1:N)]
	for iii in elemSequence
		for vk in g.adjacency_list[iii]
			if g.ordering[vk] < g.ordering[iii]
				r = vk
				while (ancestor[r] != 0) && (ancestor[r] != iii)
					t = ancestor[r]
					ancestor[r] = iii
					r = t
				end
				if ancestor[r] == 0
					ancestor[r] = iii
					par[r] = iii
				end
			end
		end
	end

	return par
end

# simplified version of my own elimination tree algorithm with simplified data structure (fastest)
function etree(g::Graph)
	N = number_of_vertices(g)
	par = zeros(Int64, N)
	# loop over Vertices of graph
	for i=1:N
		value = i
		# number of i-neighbors with order higher than order of node i
		par_ = find_parentDirect(g, i)
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
		# check degree condition for all of v's childs
		for (iii, w) in enumerate(child[v])
			# if not deg+(v) > deg+(w) - 1 for a certain w, set u to be w in snd(u), add v to snd(u)
			if !check_degree_condition(v, w, degrees)
				sn_ind[w] < 0 ? (u = w) : (u = sn_ind[w])
				sn_ind[v] = u
				child_ind = iii
				break
			end
		end

		# if v is still a rep vertex (i.e. above loop didnt findall a child that fulfilled degree condition)
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

function find_cliques(g::Graph,snodes::Array{Int64,1}, snptr::Array{Int64,1}, supernode_par::Array{Int64,1}, post::Array{Int64,1})
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
		adjPlus = find_higher_neighbors_sorted(g, vRep, postInv)
		deg = length(adjPlus) + 1
		cliques = [cliques; vRep; adjPlus]
		chptr[iii] = jjj
		nBlk[iii] = Base.power_by_squaring(length(adjPlus) + 1, 2)
		jjj += deg
	end

	return cliques, chptr, nBlk

end

function find_parentDirect(g::Graph, v::Int64)
	order = g.ordering[v]
	neighbors = g.adjacency_list[v]
	higherOrders = filter(x -> x > order, g.ordering[neighbors])
	if length(higherOrders) > 0
		return g.reverse_ordering[minimum(higherOrders)]
	else
		return 0
	end
end




