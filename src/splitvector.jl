using Base: @propagate_inbounds

const SplitView{T} = SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}

struct SplitVector{T <: AbstractFloat} <: AbstractVector{T}

	#contiguous array of source data
	data::Vector{T}
	#array of data views of type Vector{T}, with indexing of
	#each view assumed to be over a contiguous range in the typedef
	views::Vector{SplitView{T}}
	split_by::AbstractConvexSet{T}   #sets that generated the views

	#constructor (composite set)
	function SplitVector{T}(x::Vector{T}, C::CompositeConvexSet{T}) where{T}

		#check for compatibility of vector to be split
		#and sum of all of the cone sizes
		@assert C.dim == length(x)

		views = Vector{SplitView{T}}(undef, num_subsets(C))

		# loop over the sets and create views
		sidx = 0
		for i = eachindex(C.sets)
			rng = (sidx + 1):(sidx + C.sets[i].dim)
			views[i] = view(x, rng)
			sidx += C.sets[i].dim
		end
		return new(x, views, C)
	end
	#constructor (non-composite)
	function SplitVector{T}(x::Vector{T}, C::AbstractConvexSet{T}) where{T}
		views = [view(x, 1:length(x))]
		return new(x, views, C)
	end
end

SplitVector(x::Vector{T}, C) where{T} = SplitVector{T}(x, C)

Base.size(A::SplitVector) = size(A.data)
Base.length(A::SplitVector) = length(A.data)
Base.IndexStyle(::Type{<:SplitVector}) = IndexLinear()
@propagate_inbounds Base.getindex(A::SplitVector, idx::Int) = getindex(A.data, idx)
@propagate_inbounds Base.setindex!(A::SplitVector, val, idx::Int) = setindex!(A.data, val, idx)

Base.iterate(A::SplitVector) = iterate(A.data)
Base.iterate(A::SplitVector,state) = iterate(A.data, state)
Base.firstindex(A::SplitVector) = 1
Base.lastindex(A::SplitVector)  = length(A.data)

Base.showarg(io::IO, A::SplitVector, toplevel) = print(io, typeof(A))
