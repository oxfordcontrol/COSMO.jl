# -------------------------------------
# abstract type defs
# -------------------------------------
abstract type AbstractConvexSet{T} end
abstract type AbstractConvexCone{T} <: AbstractConvexSet{T} end


## -------------------------------------
# Composite Convex Set type.
# -------------------------------------

#This must be forward declared since it is used in SplitVector,
#which is itself used in most of the AbstractConvexSet methods
# See https://github.com/JuliaLang/julia/issues/269

struct CompositeConvexSet{T} <:AbstractConvexSet{T}
    dim::Int
    sets::Vector{AbstractConvexSet{T}}
    function CompositeConvexSet{T}(sets::Vector{<:AbstractConvexSet{T}}) where{T}
        # do not allow nesting of composite sets
        if any(set->isa(set,CompositeConvexSet),sets)
            throw("Nesting of CompositeConvexSets not supported")
        end
        dim = sum(x->x.dim,sets)
        new(dim,sets)
    end
end

include("./splitvector.jl")
include("./convexset.jl")
