
mutable struct Graph
	adjacency_list::Array{Array{Int64,1}}
	ordering::Array{Int64} # σ(v) = i
	reverse_ordering::Array{Int64} #σ^(-1)(i)

	function Graph(aL, o, ro)
		return new(aL, o, ro)
	end

	# constructor for list of zero or nonzero rows of vectorized matrix
	function Graph(rows::Array{Int64,1}, N::Int64, C::Union{PsdCone{Float64}, PsdConeTriangle{Float64}})
		# determine number of vertices of graph N
		ordering = collect(1:1:N)
		adjacency_list = [Int64[] for i = 1:N]

		row_val, col_val = row_ind_to_matrix_indices(rows, N, C)
		F = QDLDL.qdldl(sparse(row_val, col_val, ones(length(row_val))), logical = true)
	 	ordering = F.p
		#ordering = collect(1:N)
		lower_triangle_to_adjacency_list!(adjacency_list, F.L)

		reverse_ordering = zeros(size(ordering, 1))
		for i = 1:N
			reverse_ordering[Int64(ordering[i])] = i
		end
		return new(adjacency_list, ordering, reverse_ordering)
	end
end

# -------------------------------------
# FUNCTION DEFINITIONS
# -------------------------------------

# function to compare two graphs based on their adjacency_list
function equal_graphs(g1, g2)
	for iii = 1:length(g1.adjacency_list)
		if g1.adjacency_list[iii] != g2.adjacency_list[iii]
			return false
		end
	end
	return true
end

# deepcopy function for Graph struct
function Base.deepcopy(g::Graph)
	return Graph(deepcopy(g.adjacency_list), deepcopy(g.ordering), deepcopy(g.reverse_ordering))
end

# Redefinition of the show function that fires when the object is called
function Base.show(io::IO, obj::Graph)
	println(io, "\nGraph:\nAdjacency List:")

	# only print out information if Graph is small enough
	if length(obj.adjacency_list) <= 15
		for i= 1 :size(obj.adjacency_list, 1)
			println(io, "Vertex $(i): $(obj.adjacency_list[i])")
		end
		println(io, "\nOrdering σ(v) = i: $(obj.ordering)")
		println(io, "Reverse Ordering σ^-1(i) = v: $(obj.reverse_ordering)\n")
	end
	println(io, "\nNumber of vertices: $(length(obj.adjacency_list))\nNumber of edges: $(sum(map(x->length(x), obj.adjacency_list))/2)")
end


function number_of_vertices(g::Graph)
	return size(g.ordering, 1)
end

# returns the neighbor with the lowest order higher than the nodes order
function find_parent(g::Graph, higher_neighbors::Array{Int64})
	if size(higher_neighbors, 1) > 0
		return higher_neighbors[indmin(g.ordering[higher_neighbors])]
	else
		return 0
	end
end

# findall the cardinality of adj+(v) for all v in V
function higher_degrees(g::Graph)
	N = length(g.adjacency_list)
	degrees = zeros(Int64, N)
	for iii = 1:N
		order = g.ordering[iii]
		degrees[iii] = length(filter(x-> g.ordering[x] > order, g.adjacency_list[iii]))
	end
	return degrees
end

function find_higher_neighbors(g::Graph, node_number::Int64)
	order = g.ordering[node_number]
	neighbors = g.adjacency_list[node_number]
	higher_neighbors = neighbors[findall(f->f>order,g.ordering[neighbors])]
	return higher_neighbors
end

function find_higher_neighbors_sorted(g::Graph, node_number::Int64, ordering::Array{Int64,1})
	order = ordering[node_number]
	neighbors = g.adjacency_list[node_number]
	higher_neighbors = neighbors[findall(f -> f > order, ordering[neighbors])]
	sort!(higher_neighbors, by = x -> ordering[x])
	return higher_neighbors
end

# returns lists of Vertices that form the unconnected subgraphs (breath-first-search style)
function get_connected_parts(g::Graph)
	N = number_of_vertices(g)
	subgraphs = []
	visited = zeros(N)
	all_visited = false

	while !all_visited
		frontier = [findfirst(x -> x == 0, visited)]
		visited_nodes = [frontier[1]]
		visited[frontier[1]] = 1
		while size(frontier, 1) > 0
			next_frontier = []
			for u in frontier
				for v in g.adjacency_list[u]
					if visited[v] == 0
						push!(visited_nodes, v)
						visited[v] = 1
						push!(next_frontier, v)
					end
				end
			end
			frontier = next_frontier
		end
		# add Vertices of subgraph to array
		push!(subgraphs, visited_nodes)

		# if all Vertices are processed break
		if !in(0, visited)
			all_visited = true
		end
	end
	return subgraphs
end

# connects an unconnected graph by adding edges
function connect_graph!(g::Graph)
	subgraphs = get_connected_parts(g)

	# if more than one subgraph are found, add one edge between the first node of each subgraph
	if size(subgraphs, 1) > 1
		for i = 1:size(subgraphs, 1) -1
			node_subgraphA = subgraphs[i][1]
			node_subgraphB = subgraphs[i + 1][1]
			push!(g.adjacency_list[node_subgraphA], node_subgraphB)
			push!(g.adjacency_list[node_subgraphB], node_subgraphA)
		end
	end
	return nothing

end

# check if a graph is connected
function is_connected(g::Graph)
	return size(get_connected_parts(g),1) == 1
end

# check if the ordering of the graph is a perfect elimination ordering (i.e. for every v, are all higher neighbors adjacent?)
# start with lowest-order vertex v, findall lowest neighbor u of v with higher order. Then verify that w is adjacent to all higher order neighbors of v
# Algorithm has running time O(m+n)
function is_perfect_ordering(g::Graph)
	(length(g.reverse_ordering) == 0 || length(g.reverse_ordering) == 0) && error("Please provide graph with order and reverse order.")
	for v in g.reverse_ordering
		higher_neighbors = findhigher_neighbors(g, v)
		if size(higher_neighbors, 1) > 1
			u = higher_neighbors[indmin(g.ordering[higher_neighbors])]
			for w in higher_neighbors
				if w != u && !any(x ->x == u, g.adjacency_list[w])
					return false
				end
			end
		end
	end
	return true
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

function lower_triangle_to_adjacency_list!(alist::Array{Array{Int64, 1}, 1}, L::SparseMatrixCSC{Float64, Int64})
	N = length(alist)
	j = 1
	for (ind, r) in enumerate(L.rowval)
		i = r
		if j != N && !(ind < L.colptr[j + 1])
			j += 1
		end
		push!(alist[i], j)
		push!(alist[j], i)
	end
end


