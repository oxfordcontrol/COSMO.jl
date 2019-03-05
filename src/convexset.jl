import Base: showarg, eltype


# ----------------------------------------------------
# Zero cone
# ----------------------------------------------------
"""
    ZeroSet(dim)

Creates the zero set ``\\{ 0 \\}^{dim}`` of dimension `dim`. If `x` ∈ `ZeroSet` then all entries of x are zero.
"""
struct ZeroSet{T} <: AbstractConvexCone{T}
	dim::Int
	function ZeroSet{T}(dim::Int) where {T}
		dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
	end
end
ZeroSet(dim) = ZeroSet{DefaultFloat}(dim)


function project!(x::AbstractVector{T}, ::ZeroSet{T}) where{T}
	x .= zero(T)
	return nothing
end

function in_dual(x::AbstractVector{T}, ::ZeroSet{T}, tol::T) where{T}
	true
end

function in_pol_recc(x::AbstractVector{T}, ::ZeroSet{T}, tol::T) where{T}
	!any( x-> (abs(x) > tol), x)
end

function scale!(::ZeroSet{T}, ::SplitView{T}) where{T}
	return nothing
end

function rectify_scaling!(E,work,set::ZeroSet{T}) where{T}
	return false
end

# ----------------------------------------------------
# Nonnegative orthant
# ----------------------------------------------------
"""
    Nonnegatives(dim)

Creates the nonnegative orthant ``\\{ x \\in \\mathbb{R}^{dim} : x \\ge 0 \\}``  of dimension `dim`.
"""
struct Nonnegatives{T} <: AbstractConvexCone{T}
	dim::Int
	function Nonnegatives{T}(dim::Int) where {T}
		dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
	end
end
Nonnegatives(dim) = Nonnegatives{DefaultFloat}(dim)

function project!(x::AbstractVector{T}, C::Nonnegatives{T}) where{T}
	x .= max.(x, zero(T))
	return nothing
end

function in_dual(x::AbstractVector{T}, ::Nonnegatives{T}, tol::T) where{T}
	!any( x-> (x < -tol), x)
end

function in_pol_recc(x::AbstractVector{T}, ::Nonnegatives{T}, tol::T) where{T}
	!any( x-> (x > tol), x)
end

function scale!(cone::Nonnegatives{T}, ::AbstractVector{T}) where{T}
	return nothing
end

function rectify_scaling!(E, work, set::Nonnegatives{T}) where{T}
	return false
end

# ----------------------------------------------------
# Second Order Cone
# ----------------------------------------------------
"""
    SecondOrderCone(dim)

Creates the second-order cone (or Lorenz cone) ``\\{ (t,x) \\in \\mathrm{R}^{dim} : || x ||_2  \\leq t \\}``.
"""
struct SecondOrderCone{T} <: AbstractConvexCone{T}
	dim::Int
	function SecondOrderCone{T}(dim::Int) where {T}
		dim >= 0 ? new(dim) : throw(DomainError(dim, "dimension must be nonnegative"))
	end
end
SecondOrderCone(dim) = SecondOrderCone{DefaultFloat}(dim)

function project!(x::AbstractVector{T}, ::SecondOrderCone{T}) where{T}
	t = x[1]
	xt = view(x, 2:length(x))
	norm_x = norm(xt, 2)
	if norm_x <= t
		nothing
	elseif norm_x <= -t
		x[:] .= zero(T)
	else
		x[1] = (norm_x + t) / 2
		#x(2:end) assigned via view
		@. xt = (norm_x + t) / (2 * norm_x) * xt
	end
	return nothing
end

function in_dual(x::AbstractVector{T}, ::SecondOrderCone{T}, tol::T) where{T}
	@views norm(x[2:end]) <= (tol + x[1]) #self dual
end

function in_pol_recc(x::AbstractVector{T}, ::SecondOrderCone, tol::T) where{T}
	@views norm(x[2:end]) <= (tol - x[1]) #self dual
end

function scale!(cone::SecondOrderCone{T}, ::AbstractVector{T}) where{T}
	return nothing
end

function rectify_scaling!(E, work, set::SecondOrderCone{T}) where{T}
	return rectify_scalar_scaling!(E, work)
end

# ----------------------------------------------------
# Positive Semidefinite Cone
# ----------------------------------------------------
"""
    PsdCone(dim)

Creates the cone of symmetric positive semidefinite matrices ``\\mathcal{S}_+^{dim}``. The entries of the matrix `X` are stored column-by-column in the vector `x` of dimension `dim`.
Accordingly  ``X \\in \\mathbb{S}_+ \\Rightarrow x \\in \\mathcal{S}_+^{dim}``, where ``X = \\text{mat}(x)``.
"""
struct PsdCone{T} <: AbstractConvexCone{T}
	dim::Int
	sqrt_dim::Int
	function PsdCone{T}(dim::Int) where{T}
		dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
		iroot = isqrt(dim)
		iroot^2 == dim || throw(DomainError(dim, "dimension must be a square"))
		new(dim, iroot)
	end
end
PsdCone(dim) = PsdCone{DefaultFloat}(dim)

function project!(x::AbstractVector{T}, cone::PsdCone{T}) where{T}
	n = cone.sqrt_dim

    # handle 1D case
    if length(x) == 1
        x = max.(x, zero(T))
    else
        # symmetrized square view of x
        X    = reshape(x, n, n)
        symmetrize!(X)
        _project!(X)
    end
    return nothing
end

function _project!(X::AbstractMatrix)
	 # below LAPACK function does the following: s, U  = eigen!(Symmetric(X))
     s, U = LAPACK.syevr!('V', 'A', 'U', X, 0.0, 0.0, 0, 0, -1.0);
     # below BLAS function does the following: X .= U*Diagonal(max.(s, 0.0))*U'
     BLAS.gemm!('N', 'T', 1.0, U*Diagonal(max.(s, 0.0)), U, 0.0, X)
end

function in_dual(x::AbstractVector{T}, cone::PsdCone{T}, tol::T) where{T}
	n = cone.sqrt_dim
	X = reshape(x, n, n)
  return is_pos_sem_def(X, tol)
end

function in_pol_recc(x::AbstractVector{T}, cone::PsdCone{T}, tol::T) where{T}
	n = cone.sqrt_dim
	X = reshape(x, n, n)
	return is_neg_sem_def(X, tol)
end

function scale!(cone::PsdCone{T}, ::AbstractVector{T}) where{T}
	return nothing
end

function rectify_scaling!(E, work, set::PsdCone{T}) where{T}
	return rectify_scalar_scaling!(E, work)
end

# ----------------------------------------------------
# Positive Semidefinite Cone (Triangle)
# ----------------------------------------------------

# Psd cone given by upper-triangular entries of matrix
"""
    PsdConeTriangle(dim)

Creates the cone of symmetric positive semidefinite matrices. The entries of the upper-triangular part of matrix `X` are stored in the vector `x` of dimension `dim`.
A ``r \\times r`` matrix has ``r(r+1)/2`` upper triangular elements and results in a vector of ``\\mathrm{dim} = r(r+1)/2``.


### Examples
The matrix
```math
\\begin{bmatrix} x_1 & x_2 & x_4\\\\ x_2 & x_3 & x_5\\\\ x_4 & x_5 & x_6 \\end{bmatrix}
```
is transformed to the vector ``[x_1, x_2, x_3, x_4, x_5, x_6]^\\top `` with corresponding constraint  `PsdConeTriangle(6)`.

"""
struct PsdConeTriangle{T} <: AbstractConvexCone{T}
    dim::Int #dimension of vector
    sqrt_dim::Int # side length of matrix
    X::Array{T,2}

    function PsdConeTriangle{T}(dim::Int) where{T}
        dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
        side_dimension = Int(sqrt(0.25 + 2 * dim) - 0.5);

        new(dim, side_dimension, zeros(side_dimension, side_dimension))
    end
end
PsdConeTriangle(dim) = PsdConeTriangle{DefaultFloat}(dim)


function project!(x::AbstractArray, cone::PsdConeTriangle{T}) where{T}
    # handle 1D case
    if length(x) == 1
        x = max.(x,zero(T))
    else
        populate_upper_triangle!(cone.X, x, 1 / sqrt(2))
        _project!(cone.X)
        extract_upper_triangle!(cone.X, x, sqrt(2) )
    end
    return nothing
end

function in_dual(x::AbstractVector{T}, cone::PsdConeTriangle{T}, tol::T) where{T}
    n = cone.sqrt_dim
    populate_upper_triangle!(cone.X, x, 1 / sqrt(2))
    return is_pos_sem_def(cone.X, tol)
end

function in_pol_recc(x::AbstractVector{T}, cone::PsdConeTriangle{T}, tol::T) where{T}
    n = cone.sqrt_dim
    populate_upper_triangle!(cone.X, x, 1 / sqrt(2))
    Xs = Symmetric(cone.X)
    return is_neg_sem_def(cone.X, tol)
end

function scale!(cone::PsdConeTriangle{T}, ::AbstractVector{T}) where{T}
    return nothing
end

function rectify_scaling!(E, work, set::PsdConeTriangle{T}) where{T}
    return rectify_scalar_scaling!(E,work)
end

function populate_upper_triangle!(A::AbstractMatrix, x::AbstractVector, scaling_factor::Float64)
 	k = 0
  	for j in 1:size(A, 2)
     	for i in 1:j
        	k += 1
        	if i != j
        		A[i, j] = scaling_factor * x[k]
        	else
        		A[i, j] = x[k]
        	end
      	end
  	end
  	nothing
end

function extract_upper_triangle!(A::AbstractMatrix, x::AbstractVector, scaling_factor::Float64)
	k = 0
  	for j in 1:size(A, 2)
     	for i in 1:j
        	k += 1
        	if i != j
        		x[k] = scaling_factor * A[i, j]
        	else
        		x[k] = A[i, j]
        	end
      	end
  	end
	nothing
end

# ----------------------------------------------------
# Box
# ----------------------------------------------------
struct Box{T} <: AbstractConvexSet{T}
	dim::Int
	l::Vector{T}
	u::Vector{T}
	function Box{T}(dim::Int) where{T}
		dim >= 0 || throw(DomainError(dim, "dimension must be nonnegative"))
		l = fill!(Vector{T}(undef, dim), -Inf)
		u = fill!(Vector{T}(undef, dim), +Inf)
		new(dim, l, u)
	end
	function Box{T}(l::Vector{T}, u::Vector{T}) where{T}
		length(l) == length(u) || throw(DimensionMismatch("bounds must be same length"))
        _box_check_bounds(l,u)
        #enforce consistent bounds
		new(length(l), l, u)
	end
end
Box(dim) = Box{DefaultFloat}(dim)
Box(l, u) = Box{DefaultFloat}(l, u)

function _box_check_bounds(l,u)
    for i in eachindex(l)
        l[i] > u[i] && error("Box set: inconsistent lower/upper bounds specified at index i = ", i, ": l[i] = ",l[i],", u[i] = ",u[i])
    end
end

function project!(x::AbstractVector{T}, box::Box{T}) where{T}
	@. x = clip(x, box.l, box.u)
	return nothing
end


function support_function(x::AbstractVector{T}, B::Box{T}, tol::T) where{T}
    s = 0.
    for i in eachindex(x)
        s+= ( abs(x[i] > tol) && x[i] > 0) ? x[i]*B.u[i] : x[i]*B.l[i]
    end
    return s
end


function in_pol_recc(x::AbstractVector{T}, B::Box{T}, tol::T) where{T}
    !any(XU -> (XU[2] == Inf && XU[1] > tol), zip(x,B.u)) && !any(XL -> (XL[2] == -Inf && XL[1] < -tol), zip(x,B.l))
end

function scale!(box::Box{T}, e::AbstractVector{T}) where{T}
	@. box.l = box.l * e
	@. box.u = box.u * e
	return nothing
end

function rectify_scaling!(E, work, box::Box{T}) where{T}
	return false #no correction needed since we can scale componentwise
end


# ----------------------------------------------------
# Composite Set
# ----------------------------------------------------

#struct definition is provided in projections.jl, since it
#must be available to SplitVector, which in turn must be
#available for most of the methods here.

CompositeConvexSet(args...) = CompositeConvexSet{DefaultFloat}(args...)

function project!(x::SplitVector{T}, C::CompositeConvexSet{T}) where{T}
	@assert x.split_by === C
	foreach(xC -> project!(xC[1], xC[2]), zip(x.views, C.sets))
	return nothing
end

function support_function(x::SplitVector{T}, C::CompositeConvexSet{T}, tol::T) where{T}
	sum(xC -> support_function(xC[1], xC[2], tol),zip(x.views, C.sets))
end

function in_pol_recc(x::SplitVector{T}, C::CompositeConvexSet{T}, tol::T) where{T}
	all(xC -> in_pol_recc(xC[1], xC[2], tol), zip(x.views, C.sets))
end

function scale!(C::CompositeConvexSet{T}, e::SplitVector{T}) where{T}
	@assert e.split_by === C
	for i = eachindex(C.sets)
		scale!(C.sets[i], e.views[i])
	end
end

function rectify_scaling!(E::SplitVector{T},
	work::SplitVector{T},
	C::CompositeConvexSet{T}) where {T}
	@assert E.split_by === C
	@assert work.split_by === C
	any_changed = false
	for i = eachindex(C.sets)
		any_changed |= rectify_scaling!(E.views[i], work.views[i], C.sets[i])
	end
	return any_changed
end

#-------------------------
# general AbstractConvexCone operations
#-------------------------

# sup_{z in K_tilde_b = {-K} x {b} } <z,δy> = { <y,b> ,if y in Ktilde_polar
#                                                 +∞   ,else}

function support_function(y::SplitView{T}, cone::AbstractConvexCone{T}, tol::T) where{T}
	in_dual(-y, cone, tol) ? 0. : Inf;
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

function get_subset(C::AbstractConvexSet, idx::Int)
	idx == 1 || throw(DimensionMismatch("Input only has 1 subset (itself)"))
	return C
end
get_subset(C::CompositeConvexSet, idx::Int) = C.sets[idx]

function rectify_scalar_scaling!(E, work)
	tmp = mean(E)
	work .= tmp ./ E
	return true
end

# computes the row indices of A,b for each convex set
function get_set_indices(sets::Array{COSMO.AbstractConvexSet, 1})
	sidx = 0
	indices = Array{UnitRange{Int64}, 1}(undef, length(sets))
	for i = eachindex(sets)
		indices[i] = (sidx + 1) : (sidx + sets[i].dim)
		sidx += sets[i].dim
	end
	return indices
end
