import Base: showarg, eltype

# ----------------------------------------------------
# Zero cone
# ----------------------------------------------------
struct ZeroSet{T} <: AbstractConvexCone{T}
    dim::Int
    function ZeroSet{T}(dim::Int) where {T}
        dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
    end
end
ZeroSet(dim) = ZeroSet{DefaultFloat}(dim)


function project!(x::SplitView{T},::ZeroSet{T}) where{T}
    x .= zero(T)
    return nothing
end

function indual(x::SplitView{T},::ZeroSet{T},tol::T) where{T}
    true
end

function inrecc(x::SplitView{T},::ZeroSet{T},tol::T) where{T}
    !any(x->(abs(x) > tol),x)
end

function scale!(::ZeroSet{T},::SplitView{T}) where{T}
    return nothing
end

function rectify_scaling!(E,work,set::ZeroSet{T}) where{T}
    return false
end


# ----------------------------------------------------
# Nonnegative orthant
# ----------------------------------------------------
struct Nonnegatives{T} <: AbstractConvexCone{T}
    dim::Int
    function Nonnegatives{T}(dim::Int) where {T}
        dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
    end
end
Nonnegatives(dim) = Nonnegatives{DefaultFloat}(dim)

function project!(x::SplitView{T},C::Nonnegatives{T}) where{T}
    x .= max.(x,zero(T))
    return nothing
end

function indual(x::SplitView{T},::Nonnegatives{T},tol::T) where{T}
    !any(x->(x < -tol),x)
end

function inrecc(x::SplitView{T},::Nonnegatives{T},tol::T) where{T}
    !any(x->(x > tol),x)
end

function scale!(cone::Nonnegatives{T},::SplitView{T}) where{T}
    return nothing
end

function rectify_scaling!(E,work,set::Nonnegatives{T}) where{T}
    return false
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

function project!(x::SplitView{T},::SecondOrderCone{T}) where{T}
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
    return nothing
end

function indual(x::SplitView{T},::SecondOrderCone{T},tol::T) where{T}
    @views norm(x[2:end]) <= (tol + x[1]) #self dual
end

function inrecc(x::SplitView{T},::SecondOrderCone,tol::T) where{T}
    @views norm(x[2:end]) <= (tol - x[1]) #self dual
end

function scale!(cone::SecondOrderCone{T},::SplitView{T}) where{T}
    return nothing
end

function rectify_scaling!(E,work,set::SecondOrderCone{T}) where{T}
    return rectify_scalar_scaling!(E,work)
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

function project!(x::SplitView{T},cone::PsdCone{T}) where{T}

    n = cone.sqrtdim

    # handle 1D case
    if length(x) == 1
        x = max.(x,zero(T))
    else
        # symmetrized square view of x
        X    = reshape(x,n,n)
        X[:] = 0.5*(X+X')
        # compute eigenvalue decomposition
        # then round eigs up and rebuild
        s,U  = eigen!(X)
        floorsqrt!(s,0.)
        rmul!(U,Diagonal(s))
        mul!(X, U, U')
    end
    return nothing
end

function indual(x::SplitView{T},cone::PsdCone{T},tol::T) where{T}
    n = cone.sqrtdim
    X = reshape(x,n,n)
    return ( minimum(real(eigvals(X))) >= -tol )
end

function inrecc(x::SplitView{T},cone::PsdCone{T},tol::T) where{T}
    n = cone.sqrtdim
    X = reshape(x,n,n)
    return ( maximum(real(eigvals(X))) <= +tol )
end

function scale!(cone::PsdCone{T},::SplitView{T}) where{T}
    return nothing
end

function rectify_scaling!(E,work,set::PsdCone{T}) where{T}
    return rectify_scalar_scaling!(E,work)
end

function floorsqrt!(s::Array,floor::Real)
    @.s  = sqrt(max(floor,s))
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
        l = fill!(Vector{T}(undef,dim),-Inf)
        u = fill!(Vector{T}(undef,dim),+Inf)
        new(dim,l,u)
    end
    function Box{T}(l::Vector{T},u::Vector{T}) where{T}
        length(l) == length(u) || throw(DimensionMismatch("bounds must be same length"))
        new(length(l),l,u)
    end
end
Box(dim) = Box{DefaultFloat}(dim)
Box(l,u) = Box{DefaultFloat}(l,u)

function project!(x::SplitView{T},box::Box{T}) where{T}
    @. x = clip(x,box.l,box.u)
    return nothing
end

function indual(x::SplitView{T},box::Box{T},tol::T) where{T}
    l = box.l
    u = box.u
    for i in eachindex(x)
        if x[i] >= l[i]-tol || x[i] <= u[i]+tol
            return false
        end
    end
    return true
end

function inrecc(x::SplitView{T},::Box{T},tol::T) where{T}
    true
end

function scale!(box::Box{T},e::SplitView{T}) where{T}
    @. box.l = box.l * e
    @. box.u = box.u * e
    return nothing
end

function rectify_scaling!(E,work,box::Box{T}) where{T}
    return false #no correction needed
end


# ----------------------------------------------------
# Composite Set
# ----------------------------------------------------

#struct definition is provided in projections.jl, since it
#must be available to SplitVector, which in turn must be
#available for most of the methods here.

CompositeConvexSet(args...) = CompositeConvexSet{DefaultFloat}(args...)

function project!(x::SplitVector{T},C::CompositeConvexSet{T}) where{T}
    @assert x.splitby === C
    foreach(xC->project!(xC[1],xC[2]),zip(x.views,C.sets))
    return nothing
end

function indual(x::SplitVector{T},C::CompositeConvexSet{T},tol::T) where{T}
    all(xC -> indual(xC[1],xC[2],tol),zip(x.views,C.sets))
end

function inrecc(x::SplitVector{T},C::CompositeConvexSet{T},tol::T) where{T}
    all(xC -> inrecc(xC[1],xC[2],tol),zip(x.views,C.sets))
end

function scale!(C::CompositeConvexSet{T},e::SplitVector{T}) where{T}
    @assert e.splitby === C
    for i = eachindex(C.sets)
        scale!(C.sets[i],e.views[i])
    end
end

function rectify_scaling!(E::SplitVector{T},
                          work::SplitVector{T},
                          C::CompositeConvexSet{T}) where {T}
    @assert E.splitby === C
    @assert work.splitby === C
    any_changed = false
    for i = eachindex(C.sets)
        any_changed |= rectify_scaling!(E.views[i],work.views[i],C.sets[i])
    end
    return any_changed
end

#-------------------------
# generic set operations
#-------------------------
# function Base.showarg(io::IO, C::AbstractConvexSet{T}, toplevel) where{T}
#    print(io, typeof(C), " in dimension '", A.dim, "'")
# end

eltype(::AbstractConvexSet{T}) where{T} = T
num_subsets(C::AbstractConvexSet{T}) where{T}  = 1
num_subsets(C::CompositeConvexSet{T}) where{T} = length(C.sets)

function getsubset(C::AbstractConvexSet,idx::Int)
    idx == 1 || throw(DimensionMismatch("Input only has 1 subset (itself)"))
    return C
end
getsubset(C::CompositeConvexSet,idx::Int) = C.sets[idx]

function rectify_scalar_scaling!(E,work)
    tmp = mean(E)
    work .= tmp./E
    return true
end
