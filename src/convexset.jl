import Base: showarg, eltype

# -------------------------------------
# abstract type defs
# -------------------------------------
abstract type AbstractConvexSet{T<:AbstractFloat} end
abstract type AbstractConvexCone{T} <: AbstractConvexSet{T} end
const DefaultFloat = Float64


# ----------------------------------------------------
# Zero cone
# ----------------------------------------------------
struct ZeroCone{T} <: AbstractConvexCone{T}
    dim::Int
    function ZeroCone{T}(dim::Int) where {T}
        dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
    end
end
ZeroCone(dim) = ZeroCone{DefaultFloat}(dim)


function project!(x::AbstractVector{T},::ZeroCone{T}) where{T}
    x .= zero(T)
end

function indual(x::AbstractVector{T},::ZeroCone{T},tol::T) where{T}
    true
end

function inrecc(x::AbstractVector{T},::ZeroCone{T},tol::T) where{T}
    !any(x->(abs(x) > tol),x)
end

function scale!(::ZeroCone{T}) where{T}
    false
end


# ----------------------------------------------------
# Nonnegative orthant
# ----------------------------------------------------
struct NonnegativeOrthant{T} <: AbstractConvexCone{T}
    dim::Int
    function NonnegativeOrthant{T}(dim::Int) where {T}
        dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
    end
end
NonnegativeOrthant(dim) = NonnegativeOrthant{DefaultFloat}(dim)

function project!(x::AbstractVector{T},::NonnegativeOrthant{T}) where{T}
    @. x = max(x,zero(T))
end

function indual(x::AbstractVector{T},::NonnegativeOrthant{T},tol::T) where{T}
    !any(x->(x < -tol),x)
end

function inrecc(x::AbstractVector{T},::NonnegativeOrthant{T},tol::T) where{T}
    !any(x->(x < tol),x)
end

function scale!(cone::NonnegativeOrthant{T}) where{T}
    false
end



# ----------------------------------------------------
# Second Order Cone
# ----------------------------------------------------
struct SecondOrderCone{T} <: AbstractConvexCone{T}
    dim::Int
    function SecondOrderCone{T}(dim::Int) where {T}
        dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
    end
end
SecondOrderCone(dim) = SecondOrderCone{DefaultFloat}(dim)

function project!(x::AbstractVector{T},::SecondOrderCone{T}) where{T}
    t = x[1]
    xt = view(x,2:length(x))
    normX = norm(xt,2)
    if normX <= t
        nothing
    elseif normX <= -t
        x[:] .= zero(T)
    else
        x[1] = (normX+t)/2
        #x(2:end) assigned via view
        @.xt = (normX+t)/(2*normX)*xt
    end
    return x
end

function indual(x::AbstractVector{T},::SecondOrderCone{T},tol::T) where{T}
    @views norm(x[2:end]) <= (tol + x[1]) #self dual
end

function inrecc(x::AbstractVector{T},::SecondOrderCone,tol::T) where{T}
    @views norm(x[2:end]) <= (tol - x[1]) #self dual
end

function scale!(cone::SecondOrderCone{T}) where{T}
    true
end




# ----------------------------------------------------
# Positive Semidefinite Cone
# ----------------------------------------------------
struct PsdCone{T} <: AbstractConvexCone{T}
    dim::Int
    sqrtdim::Int
    function PsdCone{T}(dim::Int) where{T}
        dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
        iroot = isqrt(dim)
        iroot^2 == dim || throw(DomainError(x, "dimension must be a square"))
        new(dim,iroot)
    end
end
PsdCone(dim) = PsdCone{DefaultFloat}(dim)

function project!(x::AbstractVector{T},cone::PsdCone{T}) where{T}

    n = cone.sqrtdim

    # handle 1D case
    if size(x,1) == 1
        x = max.(x,0.)
    else
        # symmetrized square view of x
        X    = reshape(x,n,n)
        X[:] = 0.5*(X+X')
        # compute eigenvalue decomposition
        # then round eigs up and rebuild
        s,U  = eigen!(X)
        s  = sqrt(max(zero(T),s))
        rmul!(U,Diagonal(s))
        mul!(X, U, U')
    end
    return x
end

function indual(x::AbstractVector{T},cone::PsdCone{T},tol::T) where{T}
    n = cone.sqrtdim
    X = reshape(x,n,n)
    return ( minimum(real(eigvals(X))) >= -tol )
end

function inrecc(x::AbstractVector{T},cone::PsdCone{T},tol::T) where{T}
    n = cone.sqrtdim
    X = reshape(x,n,n)
    return ( maximum(real(eigvals(X))) <= +tol )
end

function scale!(cone::PsdCone{T}) where{T}
    true
end



# ----------------------------------------------------
# Box
# ----------------------------------------------------
struct Box{T} <:AbstractConvexSet{T}
    dim::Int
    l::Vector{T}
    u::Vector{T}
    function Box{T}(dim::Int) where{T}
        dim >= 0 || throw(DomainError(dim, "dimension must be nonnegative"))
        l = fill!(Vector{T}(undef,10),-Inf)
        u = fill!(Vector{T}(undef,10),+Inf)
        new(dim,l,u)
    end
    function Box{T}(l::Vector{T},u::Vector{T}) where{T}
        length(l) == length(u) || throw(DimensionMismatch("bounds must be same length"))
        new(length(l),l,u)
    end
end
Box(dim) = Box{DefaultFloat}(dim)
Box(l,u) = Box{DefaultFloat}(l,u)

function project!(x::AbstractVector{T},box::Box{T}) where{T}
    @. x = clip(x,box.l,box.u)
end

function indual(x::AbstractVector{T},box::Box{T},tol::T) where{T}
    l = box.l
    u = box.u
    for i in eachindex(x)
        if x[i] >= l[i]-tol || x[i] <= u[i]+tol
            return false
        end
    end
    return true
end

function inrecc(x::AbstractVector{T},::Box{T},tol::T) where{T}
    true
end

function scale!(box::Box{T}) where{T}
    false,box.l,box.u
end

# ----------------------------------------------------
# Composite Set
# ----------------------------------------------------

struct CompositeConvexSet{T} <:AbstractConvexSet{T}
    dim::Int
    sets::Vector{AbstractConvexSet{T}}
    function CompositeConvexSet{T}(sets::Vector{AbstractConvexSet{T}}) where{T}

        # do not allow nesting of composite sets
        if any(set->isa(set,CompositeConvexSet),sets)
            throw("Nesting of CompositeConvexSets not supported")
        end
        dim = sum(x->x.dim,sets)
        new(dim,sets)
    end
end
CompositeConvexSet(sets::Vector{AbstractConvexSet{T}}) where{T} = CompositeConvexSet{T}(sets)

# Functions w.r.t. CompositeConvexSets should only ever
# operate on vectors of type SplitVector.  Use
# AbstractVector{T} here instead, b/c SplitVector
# is not defined yet.
# See https://github.com/JuliaLang/julia/issues/269

function project!(x::AbstractVector{T},C::CompositeConvexSet{T}) where{T}
    @assert x.projectsto == C
    for i = eachindex(C.sets)
        project!(x.views[i],C.sets[i])
    end
end

function indual(x::AbstractVector{T},C::CompositeConvexSet{T},tol::T) where{T}
    all(xC -> indual(xC[1],xC[2],tol),zip(x,C))
end

function inrecc(x::AbstractVector{T},C::CompositeConvexSet{T},tol::T) where{T}
    all(xC -> inrecc(xC[1],xC[2],tol),zip(x,C))
end

function scale!(C::CompositeConvexSet{T}) where{T}
    nothing
end

#-------------------------
# generic set operations
#-------------------------
# function Base.showarg(io::IO, C::AbstractConvexSet{T}, toplevel) where{T}
#    print(io, typeof(C), " in dimension '", A.dim, "'")
# end

Base.eltype(::AbstractConvexSet{T}) where{T} = T
num_subsets(C::AbstractConvexSet)  = 1
num_subsets(C::CompositeConvexSet) = length(C.sets)
