"""
    compute_reduced_clique_graph!(sep::Array{Set{Int64}, 1}, snd::Array{Set{Int64}, 1})

Compute the reduced clique graph (union of all clique trees) given an initial clique tree defined by its supernodes `snd` and separator `sep` sets.

We are using the algorithm described in **Michel Habib and Juraj Stacho - Polynomial-time algorithm for the leafage ofchordal graphs (2009)**, which
computes the reduced clique graph in the following way:
1. Sort all minimal separators by size
2. Initialise graph CG(R) with cliques as nodes and no edges
3. for largest unprocessed separator S and
    |  add an edge between any two cliques C1 and C2 if they both contain S and are in different connected   components   of CG(R) and store in `edges`.
    |  Compute an edge weight used for merge decision and store in `val`.
    |  Store the index of the separator which is the intersection C1 ∩ C2 in `iter`
   end
"""
function compute_reduced_clique_graph!(sep::Array{Set{Int64}, 1}, snd::Array{Set{Int64}, 1})
    # loop over separators by decreasing cardinality
    sort!(sep, by = x -> length(x), rev = true)

    # edges = Array{Tuple{Int64, Int64}, 1}()    # a list of edges in the reduced clique graph, higher clique index first
    rows = Int64[]
    cols = Int64[]
    # inter = Array{Int64, 1}()                  # the index of the separator which corresponds to the intersection of the two cliques

    for (k, separator) in enumerate(sep)
        # find cliques that contain the separator
        clique_ind = findall(x -> separator ⊆ x, snd)

        # we compute the separator graph (see Habib, Stacho - Reduced clique graphs of chordal graphs) to analyse connectivity
        # we represent the separator graph H by a hashtable
        H = separator_graph(clique_ind, separator, snd)
        # find the connected components of H
        components = find_components(H, clique_ind)
        # for each pair of cliques that contain the separator, add an edge to the reduced clique tree if they are in unconnected components
        for pair in subsets(clique_ind, Val{2}())
            if is_unconnected(pair, components)
                # push!(edges, (max(pair...), min(pair...))) #add edge
                push!(rows, max(pair...))
                push!(cols, min(pair...))
                # push!(inter, k) # store intersection
            end
        end

    end
    return rows, cols#, inter
end

"Check whether the `pair` of cliques are in different `components`."
function is_unconnected(pair::Tuple{Int64, Int64}, components::Array{Set{Int64}, 1})
    component_ind = findfirst(x -> pair[1] ∈ x, components)
    return pair[2] ∉ components[component_ind]
end

"Find the separator graph H given a separator and the relevant index-subset of cliques."
function separator_graph(clique_ind::Array{Int64,1}, separator::Set{Int64}, snd::Array{Set{Int64}, 1})

    # make the separator graph using a hash table
    # key: clique_ind --> edges to other clique indices
    H = Dict{Int64, Array{Int64, 1}}()

    for pair in subsets(clique_ind, Val{2}())
        ca = pair[1]
        cb = pair[2]
         # if intersect_dim(snd[ca], snd[cb]) > length(separator)
        if !inter_equal(snd[ca], snd[cb], separator)
            if haskey(H, ca)
                push!(H[ca], cb)
            else
                H[ca] = [cb]
            end
            if haskey(H, cb)
                push!(H[cb], ca)
            else
                H[cb] = [ca]
            end
        end
    end
    # add unconnected cliques
    for v in clique_ind
        !haskey(H, v) && (H[v] = Int64[])
    end
    return H
end

"Find connected components in undirected separator graph represented by `H`."
function find_components(H::Dict{Int64, Array{Int64, 1}}, clique_ind::Array{Int64,1})
    visited = Dict{Int64, Bool}(v => false for v in clique_ind)
    components = Array{Set{Int64}, 1}()
    for v in clique_ind
        if visited[v] == false
            component = Set{Int64}()
            push!(components, DFS_hashtable!(component, v, visited, H))
        end
    end
    return components
end

"Depth first search on a hashtable `H`."
function DFS_hashtable!(component::Set{Int64}, v::Int64, visited::Dict{Int64, Bool}, H::Dict{Int64, Array{Int64, 1}})
    visited[v] = true
    push!(component, v)
    for n in H[v]
        if visited[n] == false
            DFS_hashtable!(component, n, visited, H)
        end
    end
    return component
end


"Return the number of elements in s ∩ s2."
function intersect_dim(s::Set, s2::Set)
    dim = 0
    for e in s
        e in s2 && (dim += 1)
    end
    return dim
end

"Check if s ∩ s2 == s3."
function inter_equal(s::Set{Int64}, s2::Set{Int64}, s3::Set{Int64})
    dim = 0
    len_s3 = length(s3)
    if length(s) < length(s2)
        sa = s
        sb = s2
    else
        sa = s2
        sb = s
    end
    for e in sa
        if e in sb
            dim += 1
            dim > len_s3 && return false
            e in s3 || return false
        end
    end
    return true
end


"Given a list of edges, return an adjacency hash-table `table` with nodes from 1 to `num_vertices`."
function compute_adjacency_table(edges::SparseMatrixCSC{Float64, Int64}, num_vertices::Int64)
    table = Dict(i => Set{Int64}() for i = 1:num_vertices)
    r = edges.rowval
    c = edges.colptr
     for col = 1:num_vertices
         for i in c[col]:c[col+1]-1
             row = r[i]
             push!(table[row], col)
            push!(table[col], row)
        end
     end
     return table
end

"Check whether `edge` is permissible for a merge. An edge is permissible if for every common neighbor N, C_1 ∩ N == C_2 ∩ N or if no common neighbors exist."
function ispermissible(edge::Tuple{Int64, Int64}, adjacency_table::Dict{Int64, Set{Int64}}, snd::Array{Set{Int64}, 1})
    c_1 = edge[1]
    c_2 = edge[2]
    common_neighbors = intersect(adjacency_table[c_1], adjacency_table[c_2])
    # N.B. This can be made faster by first checking whether the sizes of the intersection are the same before allocating anything
    for neighbor in common_neighbors
        intersect(snd[c_1], snd[neighbor]) != intersect(snd[c_2], snd[neighbor]) && return false
    end
    return true
end


heapleft(i::Integer) = 2i
heapright(i::Integer) = 2i + 1
heapparent(i::Integer) = div(i, 2)


mutable struct BinaryHeapMax{T}
    order::Vector{Int64}
    weights::Vector{T}

    function BinaryHeapMax{T}(weights::Array{T}) where {T <: Real}
        order = collect(1:length(weights))
        build_max_heap!(order, weights)
        return new(order, weights) # weights is a reference on edges.nzval
    end
end
BinaryHeapMax(weights::Array{T}) where {T <: Real} = BinaryHeapMax{T}(weights)

function build_max_heap!(arr::Vector{Int64}, vals::Vector{T}) where {T}
    for i = div(length(arr), 2):-1:1
        max_heapify!(arr, vals, i)
    end
end

function max_heapify!(arr::Vector{Int64}, vals::Vector{T}, i::Int64) where {T}
    println("Start with i=$(i)")
    left = heapleft(i)
    right = heapright(i)
    largest = copy(i)
    @show(left, right, largest)
    @show(arr, vals)
    if left <= length(arr) && vals[arr[left]] > vals[arr[largest]]
        largest = left
    end
    if right <= length(arr) && vals[arr[right]] > vals[arr[largest]]
        largest = right
    end

    if largest != i
        ai = arr[i]
        arr[i] = copy(arr[largest])
        arr[largest] = ai
        @show(i, largest, arr)
        println("Next step:")
        max_heapify!(arr, vals, largest)
    else
        println("Nope")
    end
    return nothing
end
