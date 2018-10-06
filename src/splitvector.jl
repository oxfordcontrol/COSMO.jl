using Base: @propagate_inbounds

abstract type AbstractConvexSet{T<:AbstractFloat} end

struct SplitVector{T<:AbstractFloat} <: AbstractVector{T}

    #contiguous array of source data
    data::Vector{T}
    #array of data views of type Vector{T}, with indexing of
    #each view assumed to be over a contiguous range
    views::Array{SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}}
    projectsto::AbstractConvexSet{T}   #sets that generated the views

    #constructor (composite set)
    function SplitVector{T}(x::Vector{T},C::CompositeConvexSet{T}) where{T}

        #check for compatibility of vector to be split
        #and sum of all of the cone sizes
        @assert sum(set->set.dim,C.sets) == length(x)

        #I want an array of views. The actual type
        #of a view is convoluted, so just make one
        #and test it directly.  Use 'similar' in case x
        #is length zero for some reason
        vtype = typeof(view(similar(x,1),1:1))
        views = Array{vtype}(undef,num_subsets(C))

        # loop over the sets and create views
        sidx = 0
        for i = eachindex(C.sets)
            rng = (sidx+1):(sidx + C.sets[i].dim)
            views[i] = view(x,rng)
            sidx += C.sets[i].dim
        end
        return new(x,views,C)
    end
    #constructor (non-composite)
    function SplitVector{T}(x::Vector{T},C::AbstractConvexSet{T}) where{T}
        views = [view(x,:)]
        return new(x,views,C)
    end
end

SplitVector(x::Vector{T},C) where{T} = SplitVector{T}(x,C)

Base.size(A::SplitVector) = size(A.data)
Base.length(A::SplitVector) = length(A.data)
Base.IndexStyle(::Type{<:SplitVector}) = IndexLinear()
@propagate_inbounds Base.getindex(A::SplitVector, idx::Int) = getindex(A.data,idx)
@propagate_inbounds Base.setindex!(A::SplitVector, val, idx::Int) = setindex!(A.data,val,idx)

Base.iterate(A::SplitVector) = iterate(A.data)
Base.iterate(A::SplitVector,state) = iterate(A.data,state)
Base.firstindex(A::SplitVector) = 1
Base.lastindex(A::SplitVector)  = length(A.data)

Base.showarg(io::IO, A::SplitVector, toplevel) = print(io, typeof(A))
