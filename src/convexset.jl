import Base: showarg, eltype
using Arpack, LinearMaps, LinearAlgebra

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
mutable struct PsdCone{T} <: AbstractConvexCone{T}
    dim::Int
    sqrtdim::Int
    positive_subspace::Bool
    Z::AbstractMatrix{Float64}  # Ritz vectors
    λ::AbstractVector{Float64}  # Ritz values
    λ_rem::T
    z_rem::AbstractVector{Float64}
    buffer_size::Int
    iter_number::Int
    function PsdCone{T}(dim::Int) where{T}
        dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
        iroot = isqrt(dim)
        iroot^2 == dim || throw(DomainError(x, "dimension must be a square"))
        new(dim,iroot,true,zeros(T,iroot,0),zeros(T,iroot),0.0, randn(T,iroot), 3, 0)
    end
end
PsdCone(dim) = PsdCone{DefaultFloat}(dim)

function project_to_nullspace(X::AbstractArray, x::AbstractVector, tmp::AbstractVector)
    # Project x to the nullspace of X', i.e. x .= (I - X*X')*x
    # tmp is a vector used in intermmediate calculations
    T = eltype(X)
    BLAS.gemv!('T', one(T), X, x, zero(T), tmp)
    BLAS.gemv!('N', -one(T), X, tmp, one(T), x)
end

function estimate_λ_rem(X::AbstractArray, U::AbstractArray, n::Int, x0::AbstractVector)
	# Estimates largest eigenvalue of the Symmetric X on the subspace we discarded
	# Careful, we need to enforce all iterates to be orthogonal to the range of U
    T = eltype(X)
    tmp = zeros(T, size(U, 2))
    function custom_mul!(y::AbstractVector, x::AbstractVector)
        # Performs y .= (I - U*U')*X*x
        # y .= X*x - U*(U'*(X*x))
        BLAS.symv!('U', one(T), X, x, zero(T), y)
        project_to_nullspace(U, y, tmp)
	end
    project_to_nullspace(U, x0, tmp)
    A = LinearMap{T}(custom_mul!, size(X, 1); ismutating=true, issymmetric=true)
    (λ_rem, v_rem, nconv, niter, nmult, resid) = eigs(A, nev=n,
        ncv=20, ritzvec=true, which=:LR, tol=.1, v0=x0)
    return λ_rem, v_rem
end

function generate_subspace(X::AbstractArray, cone::PsdCone) 
    W = Array(qr([cone.Z X*cone.Z]).Q)
    XW = X*W
    return W, XW
end

function propagate_subspace(Z, XZ, ZXZ, λ, tol=1e-5)
    res_norms = colnorms(XZ - Z*Diagonal(λ)) 
    XV = XZ[:, res_norms .> tol]
    VXV = ZXZ[:, res_norms .> tol]
    VX2V = XV'*XV;
    return Z, XV, VXV, VX2V
end

function generate_subspace_chol(X::AbstractArray, cone::PsdCone)
    XZ = X*cone.Z; ZXZ = cone.Z'*XZ
    V, XV, VXV, VX2V = propagate_subspace(cone.Z, XZ, ZXZ, cone.λ)
    if size(VXV, 2) == 0
        return cone.Z, XZ
    end
    R = cholesky([I VXV; VXV' VX2V]; check=false)
    if !issuccess(R) || abs(logdet(R)) > 100
        #=
        λ, U = eigen(Symmetric(ZXZ))
        V, XV, VXV, VX2V = propagate_subspace(cone.Z, XZ*U, ZXZ*U, λ)
        if size(VXV, 2) == 0
            return cone.Z, XZ
        end
        R = cholesky([I VXV; VXV' VX2V])
        =#
        @show "QR!"
        W = Array(qr([cone.Z XZ]).Q)
        return W, X*W
    end
    @show logdet(R)
    W = (R.L\[V XV]')';
    @show norm(W'*W - I)
    XW = X*W

    return W, XW
end

function project!(x::AbstractArray,cone::PsdCone{T}) where{T}
    n = cone.sqrtdim
    cone.iter_number += 1

    if mod(cone.iter_number, 40) == 0 || size(cone.Z, 2) == 0 || size(cone.Z, 2) >= cone.sqrtdim/2 || n == 1
        # @show size(cone.Z)
        return project_exact!(x, cone)
    end

    if !cone.positive_subspace
        @. x = -x
    end
    X = Symmetric(reshape(x, n, n));
    X = convert(Array{Float64, 2}, X) # Convert to Float32

    W, XW = generate_subspace(X, cone)
    Xsmall = W'*XW
    l, V, first_positive, first_negative = eigen_sorted(Symmetric(Xsmall), 1e-10);
	
    # Positive Ritz pairs
    Vp = V[:, first_positive:end]
    U = W*Vp; λ = l[first_positive:end];

    buffer_idx = max(first_negative-cone.buffer_size-1,1):max(first_negative,0)
	Ub = W*V[:, buffer_idx]; λb = l[buffer_idx]

	# Projection
    Xπ = U*Diagonal(λ)*U';
    
    # Residual Calculation
    # R = XW*Vp - U*Diagonal(λ)
    R = X*U - U*Diagonal(λ)

    # Important: why not [U Ub]?
    λ_rem, z_rem = estimate_λ_rem(X, U, 2, cone.z_rem)
    λ_rem .= max.(λ_rem, 0.0)
    eig_sum = sum(λ_rem).^2 + (n - size(W, 2) - length(λ_rem))*minimum(λ_rem).^2
    @show residual = sqrt(2*norm(R, 2)^2 + eig_sum)
    #=
    @show norm(R, 2)
    @show sqrt(eig_sum)
    =#
    
    if !cone.positive_subspace
        @. Xπ -= X 
    end
    copyto!(x, Xπ)

    # It should be [U ub z_rem]
    cone.Z = [U z_rem]
    cone.z_rem = randn(size(U, 1))
    cone.λ = [λ; λ_rem]
end

function project_exact!(x::AbstractArray{T},cone::PsdCone{T}) where{T}
    n = cone.sqrtdim

    # handle 1D case
    if length(x) == 1
        x = max.(x,zero(T))
    else
        # symmetrized square view of x
        X    = reshape(x,n,n)
        @. X = (X + X')/2
        # compute eigenvalue decomposition
        # then round eigs up and rebuild
        λ, U  = eigen!(Symmetric(X))
        Up = U[:, λ .> 0]
        sqrt_λp = sqrt.(λ[λ .> 0])
        if length(sqrt_λp) > 0
            # rmul!(Up, Diagonal(sqrt_λp))
            # mul!(X, Up, Up')
            X = Up*Diagonal(λ[λ .> 0])*Up'
        else
            X .= 0
            return nothing
            #ToDo: Handle this case with lanczos
        end
        x .= reshape(X, cone.dim)
        
        # @show λ
        # Save the subspace we will be tracking
        if sum(λ .> 0) <= sum(λ .< 0)
            sorted_idx = sortperm(λ)
            cone.positive_subspace = true
            idx = findfirst(λ[sorted_idx] .> 0) # First positive index
            if isa(idx, Nothing)
                idx = length(λ) + 1
            end
        else
            sorted_idx = sortperm(-λ)
            cone.positive_subspace = false
            idx = findfirst(λ[sorted_idx] .< 0) # First positive index
            if isa(idx, Nothing)
                idx = length(λ) + 1
            end
        end
        # Take also a few vectors from the other discarted eigenspace
        idx = max(idx - cone.buffer_size, 1)
        cone.Z = U[:, sorted_idx[idx:end]]
        cone.λ = λ[sorted_idx[idx:end]]
        if sum(λ .> 0) >= sum(λ .< 0)
            cone.λ = -cone.λ
        end
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
