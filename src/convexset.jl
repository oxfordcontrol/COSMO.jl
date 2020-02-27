using UnsafeArrays
import Base: showarg, eltype
const DSYEVR_ = (BLAS.@blasfunc(dsyevr_),Base.liblapack_name)
const SSYEVR_ = (BLAS.@blasfunc(ssyevr_),Base.liblapack_name)

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


function allocate_memory!(cone::AbstractConvexSet{T}) where {T}
  return nothing
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

# ----------------------------------------------------
# Positive Semidefinite Cone
# ----------------------------------------------------

#a type to maintain internal workspace data for the BLAS syevr function
mutable struct PsdBlasWorkspace{T}
    m::Base.RefValue{BLAS.BlasInt}
    w::Vector{T}
    Z::Matrix{T}
    isuppz::Vector{BLAS.BlasInt}
    work::Vector{T}
    lwork::BLAS.BlasInt
    iwork::Vector{BLAS.BlasInt}
    liwork::BLAS.BlasInt
    info::Base.RefValue{BLAS.BlasInt}

    function PsdBlasWorkspace{T}(n::Integer) where{T}

        BlasInt = BLAS.BlasInt

        #workspace data for BLAS
        m      = Ref{BlasInt}()
        w      = Vector{T}(undef,n)
        Z      = Matrix{T}(undef,n,n)
        isuppz = Vector{BlasInt}(undef, 2*n)
        work   = Vector{T}(undef, 1)
        lwork  = BlasInt(-1)
        iwork  = Vector{BlasInt}(undef, 1)
        liwork = BlasInt(-1)
        info   = Ref{BlasInt}()

        new(m,w,Z,isuppz,work,lwork,iwork,liwork,info)
    end
end

for (syevr, elty) in
    ((DSYEVR_,:Float64),
     (SSYEVR_,:Float32))
   @eval begin
        function _syevr!(A::AbstractMatrix{$elty}, ws::PsdBlasWorkspace{$elty})

            #Float64 only support for now since we call dsyevr_ directly
            n       = size(A,1)
            ldz     = n
            lda     = stride(A,2)

                ccall($syevr, Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BLAS.BlasInt},
                Ptr{$elty}, Ref{BLAS.BlasInt}, Ref{$elty}, Ref{$elty},
                Ref{BLAS.BlasInt}, Ref{BLAS.BlasInt}, Ref{$elty}, Ptr{BLAS.BlasInt},
                Ptr{$elty}, Ptr{$elty}, Ref{BLAS.BlasInt}, Ptr{BLAS.BlasInt},
                Ptr{$elty}, Ref{BLAS.BlasInt}, Ptr{BLAS.BlasInt}, Ref{BLAS.BlasInt},
                Ptr{BLAS.BlasInt}),
                'V', 'A', 'U', n,
                A, max(1,lda), 0.0, 0.0,
                0, 0, -1.0,
                ws.m, ws.w, ws.Z, ldz, ws.isuppz,
                ws.work, ws.lwork, ws.iwork, ws.liwork,
                ws.info)
                LAPACK.chklapackerror(ws.info[])
        end
    end #@eval
end #for

function _project!(X::AbstractMatrix, ws::PsdBlasWorkspace{T}) where{T}

    #computes the upper triangular part of the projection of X onto the PSD cone

     #allocate additional workspace arrays if the ws
     #work and iwork have not yet been sized
     if ws.lwork == -1
         _syevr!(X,ws)
         ws.lwork = BLAS.BlasInt(real(ws.work[1]))
         resize!(ws.work, ws.lwork)
         ws.liwork = ws.iwork[1]
         resize!(ws.iwork, ws.liwork)
     end

	 # below LAPACK function does the following: w,Z  = eigen!(Symmetric(X))
	 	_syevr!(X, ws)
		# compute upper triangle of: X .= Z*Diagonal(max.(w, 0.0))*Z'
		rank_k_update!(X, ws)
end

function rank_k_update!(X::AbstractMatrix, ws::COSMO.PsdBlasWorkspace{T}) where{T}
  n = size(X, 1)
  X .= 0
  nnz_λ = 0
  for j = 1:length(ws.w)
    λ = ws.w[j]
    if λ > 0
      nnz_λ += 1
      @inbounds for i = 1:n
        ws.Z[i, j] = ws.Z[i, j] * sqrt(λ)
      end
    end
  end

  if nnz_λ > 0
    V = uview(ws.Z, :, (n - nnz_λ + 1):n)
    BLAS.syrk!('U', 'N', 1.0, V, 1.0, X)
  end
  return nothing
end

"""
    PsdCone(dim)

Creates the cone of symmetric positive semidefinite matrices ``\\mathcal{S}_+^{dim}``. The entries of the matrix `X` are stored column-by-column in the vector `x` of dimension `dim`.
Accordingly  ``X \\in \\mathbb{S}_+ \\Rightarrow x \\in \\mathcal{S}_+^{dim}``, where ``X = \\text{mat}(x)``.
"""
struct PsdCone{T} <: AbstractConvexCone{T}
	dim::Int
	sqrt_dim::Int
  work::PsdBlasWorkspace{T}
  tree_ind::Int64  # tree number that this cone belongs to
  clique_ind::Int64
	function PsdCone{T}(dim::Int, tree_ind::Int64, clique_ind::Int64) where{T}
		dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
		iroot = isqrt(dim)
		iroot^2 == dim || throw(DomainError(dim, "dimension must be a square"))
		new(dim, iroot, PsdBlasWorkspace{T}(iroot), tree_ind, clique_ind)
	end
end
PsdCone(dim) = PsdCone{DefaultFloat}(dim)
PsdCone{T}(dim::Int64) where{T} = PsdCone{T}(dim, 0, 0)


struct DensePsdCone{T} <: AbstractConvexCone{T}
  dim::Int
  sqrt_dim::Int
  work::PsdBlasWorkspace{T}
  function DensePsdCone{T}(dim::Int) where{T}
    dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
    iroot = isqrt(dim)
    iroot^2 == dim || throw(DomainError(dim, "dimension must be a square"))
    new(dim, iroot, PsdBlasWorkspace{T}(iroot))
  end
end
DensePsdCone(dim) = DensePsdCone{DefaultFloat}(dim)



function project!(x::AbstractVector{T}, cone::Union{PsdCone{T}, DensePsdCone{T}}) where{T}
	n = cone.sqrt_dim

    # handle 1D case
    if length(x) == 1
        x .= max(x[1], zero(T))
    else
        # symmetrized square view of x
        X    = reshape(x, n, n)
        symmetrize_upper!(X)
        _project!(X,cone.work)

        #fill in the lower triangular part
        for j=1:n, i=1:(j-1)
            X[j,i] = X[i,j]
        end
    end
    return nothing
end

# Notice that this is an in-place version that uses x as workspace
function in_dual!(x::AbstractVector{T}, cone::Union{PsdCone{T}, DensePsdCone{T}}, tol::T) where{T}
	n = cone.sqrt_dim
	X = reshape(x, n, n)
  return COSMO.is_pos_def!(X, tol)
end
in_dual(x::AbstractVector{T}, cone::Union{PsdCone{T}, DensePsdCone{T}}, tol::T) where{T} = in_dual!(copy(x), cone, tol)

function in_pol_recc!(x::AbstractVector{T}, cone::Union{PsdCone{T}, DensePsdCone{T}}, tol::T) where{T}
	n = cone.sqrt_dim
	X = reshape(x, n, n)
	return COSMO.is_neg_def!(X, tol)
end
in_pol_recc(x::AbstractVector{T}, cone::Union{PsdCone{T}, DensePsdCone{T}}, tol::T) where{T} = in_pol_recc!(copy(x), cone, tol)




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
mutable struct PsdConeTriangle{T} <: AbstractConvexCone{T}
    dim::Int #dimension of vector
    sqrt_dim::Int # side length of matrix
    X::Array{T,2}
    work::PsdBlasWorkspace{T}
    tree_ind::Int64 # tree number that this cone belongs to
    clique_ind::Int64

    function PsdConeTriangle{T}(dim::Int, tree_ind::Int64, clique_ind::Int64) where{T}
      dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
      side_dimension = Int(sqrt(0.25 + 2 * dim) - 0.5);
      new(dim, side_dimension, zeros(side_dimension, side_dimension), PsdBlasWorkspace{T}(side_dimension), tree_ind, clique_ind)

  end
end
PsdConeTriangle(dim) = PsdConeTriangle{DefaultFloat}(dim)
PsdConeTriangle{T}(dim::Int64) where{T} = PsdConeTriangle{T}(dim, 0, 0)

DecomposableCones{T} = Union{PsdCone{T}, PsdConeTriangle{T}}

mutable struct DensePsdConeTriangle{T} <: AbstractConvexCone{T}
    dim::Int #dimension of vector
    sqrt_dim::Int # side length of matrix
    X::Array{T,2}
    work::PsdBlasWorkspace{T}

    function DensePsdConeTriangle{T}(dim::Int) where{T}
        dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
        side_dimension = Int(sqrt(0.25 + 2 * dim) - 0.5);
        new(dim, side_dimension, zeros(side_dimension, side_dimension), PsdBlasWorkspace{T}(side_dimension))
    end
end
DensePsdConeTriangle(dim) = DensePsdConeTriangle{DefaultFloat}(dim)



function project!(x::AbstractArray, cone::Union{PsdConeTriangle{T}, DensePsdConeTriangle{T}}) where{T}
    # handle 1D case
    if length(x) == 1
        x .= max(x[1],zero(T))
    else
        populate_upper_triangle!(cone.X, x, 1 / sqrt(2))
        _project!(cone.X,cone.work)
        extract_upper_triangle!(cone.X, x, sqrt(2) )
    end
    return nothing
end

# Notice that we are using a (faster) in-place version that modifies the input
function in_dual!(x::AbstractVector{T}, cone::Union{PsdConeTriangle{T}, DensePsdConeTriangle{T}, PsdConeTriangleLOBPCG{T}}, tol::T) where{T}
    n = cone.sqrt_dim
    populate_upper_triangle!(cone.X, x, 1 / sqrt(2))
    return COSMO.is_pos_def!(cone.X, tol)
end
in_dual(x::AbstractVector{T}, cone::Union{PsdConeTriangle{T}, DensePsdConeTriangle{T}, PsdConeTriangleLOBPCG{T}}, tol::T) where {T} = in_dual!(x, cone, tol)

function in_pol_recc!(x::AbstractVector{T}, cone::Union{PsdConeTriangle{T}, DensePsdConeTriangle{T}, PsdConeTriangleLOBPCG{T}}, tol::T) where{T}
    n = cone.sqrt_dim
    populate_upper_triangle!(cone.X, x, 1 / sqrt(2))
    Xs = Symmetric(cone.X)
    return COSMO.is_neg_def!(cone.X, tol)
end
in_pol_recc(x::AbstractVector{T}, cone::Union{PsdConeTriangle{T}, DensePsdConeTriangle{T}, PsdConeTriangleLOBPCG{T}}, tol::T) where {T} = in_pol_recc!(x, cone, tol)


function allocate_memory!(cone::Union{PsdConeTriangle{T}, DensePsdConeTriangle{T}}) where {T}
  cone.X = zeros(cone.sqrt_dim, cone.sqrt_dim)
end

function reset_iteration_counters!(cone::AbstractConvexSet)
  nothing
end

"""
    ExponentialCone(MAX_ITERS = 100, EXP_TOL = 1e-8)

Creates the exponential cone ``\\mathcal{K}_{exp} = \\{(x, y, z) \\mid y \\geq 0 ye^{x/y} ≤ z\\} \\cup \\{ (x,y,z) \\mid   x \\leq 0, y = 0, z \\geq 0 \\}``
"""
struct ExponentialCone{T} <: AbstractConvexCone{T}
  dim::Int
  v0::Vector{T}
  MAX_ITER::Int64
  EXP_TOL::Float64

  function ExponentialCone{T}(dim = 3, MAX_ITERS = 100, EXP_TOL = 1e-8) where{T}
    new(3, zeros(T, 3), MAX_ITERS, EXP_TOL)
  end
end
ExponentialCone(args...) = ExponentialCone{DefaultFloat}(args...)


function project!(v::AbstractVector{T}, cone::ExponentialCone{T}) where{T}

  # Check the four different cases
  # 1. v in K_exp => v = v
  in_cone(v, cone, 0.) && return nothing

  # 2. -v in K_exp^* => v = 0
  if in_dual(-v, cone, 0.)
    v .= zero(T)
    return nothing
  end

  # 3. x < 0 and y < 0 => v = (x, 0, max(z, 0))
  if v[1] < 0 && v[2] < 0
    v[2] = 0.0
    v[3] = max(v[3], 0)
    return nothing
  end

  # 4. Otherwise solve the following minimisation problem
  # min_w (1/2) ||v - v0||_2^2
  # s.t.  y * exp(x/y) == z
  #       y > 0
  project_exp!(v, cone)
end

# This is a modified version of the projection code used in SCS
# https://github.com/cvxgrp/scs/blob/master/src/cones.c
# We are solving the dual problem g(λ) via a bisection method
function project_exp!(v::AbstractVector{T}, cone::ExponentialCone{T}) where{T}
  # save input vector and use v as working variable
  @. cone.v0 = v
  l, u = get_bisection_bounds(v, cone.v0, cone.EXP_TOL)

  for k = 1:cone.MAX_ITER
    λ = (u + l) / 2
    g = grad_dual!(λ, v, cone.v0, cone.EXP_TOL)
    g > 0 ? (l = λ) : (u = λ)
    u - l < cone.EXP_TOL && break
  end
end

function get_bisection_bounds(v::AbstractVector{T}, v0::Vector{T}, tol::Float64) where {T <: Real}
  l = 0.
  λ = 0.125
  g = grad_dual!(λ, v, v0, tol)
  while g > 0
    l = λ
    λ *= 2
    g = grad_dual!(λ, v, v0, tol)
  end
  u = λ
  return l, u
end

function grad_dual!(λ::T, v::AbstractVector{T}, v0::Vector{T}, tol::Float64) where {T <: Real}
  find_minimizers!(λ, v, v0, tol)
  v[2] == 0 ? (g = v[1]) : (g = v[1] + v[2] * log(v[2] / v[3]))
  return g
end

function find_minimizers!(λ::T, v::AbstractVector{T}, v0::Vector{T}, tol::Float64) where {T <: Real}
  v[3] = find_min_t(λ, v0[2], v0[3], tol)
  # s* = (t - t0) * t / λ
  v[2] = (1 / λ) * (v[3] - v0[3]) * v[3]
  # r* = r0 - λ
  v[1] = v0[1] - λ
end

# use Newton method to find minimizer t* for given λ, i.e. find the zero of
# f(t) = t * (t - t0) / lambda - s0 + λ * log( t - t0 / λ) + λ
# Define Δt = t - t0
function find_min_t(λ::T, s0::T, t0::T, tol::Float64) where {T <:Real}
  Δt = max(-t0, tol)
  for k = 1:150
    f = Δt * (Δt + t0) / λ^2 - s0 / λ + log(Δt / λ) + 1
    grad_f = (2 * Δt + t0) / λ^2 + 1 / Δt
    Δt = Δt - f / grad_f

    if (Δt <= -t0)
      Δt = -t0
      break
    elseif (Δt <= 0)
      Δt = 0
      break
    elseif abs(f) < tol
      break
    end
  end
  return Δt + t0
end

function in_cone(v::AbstractVector{T}, cone::ExponentialCone{T}, tol::T) where{T}
  x = v[1]
  y = v[2]
  z = v[3]
  return (y > 0 && y * exp(x/y) <= z + tol) || (x <= tol &&  y == 0. && z >= -tol )
end
# Kexp^* = { (x,y,z) | x < 0, -xe^(y/x) <= e^1 z } cup { (0,y,z) | y >= 0,z >= 0 }
function in_dual(v::AbstractVector{T}, cone::ExponentialCone{T}, tol::T) where{T}
  x = v[1]
  y = v[2]
  z = v[3]
  return (x < 0 && -x * exp(y / x) - exp(1) *  z <= tol) || (abs(x) <= tol && y >= -tol && z >= -tol)
end

function in_pol_recc(v::AbstractVector{T},cone::ExponentialCone{T}, tol::T) where{T}
  return in_dual(-v, cone, tol)
end


"""
    PowerCone(alpha::Float64, MAX_ITERS::Int64 = 20, POW_TOL = 1e-8)

Creates the 3-d power cone ``\\mathcal{K}_{pow} = \\{(x, y, z) \\mid x^\\alpha y^{(1-\\alpha)} \\geq  \\|z\\|, x \\geq 0, y \\geq 0 \\}`` with ``0 < \\alpha < 1``
"""
struct PowerCone{T} <: AbstractConvexCone{T}
  dim::Int
  α::Float64
  MAX_ITER::Int64
  POW_TOL::Float64

  function PowerCone{T}(alpha::Float64, MAX_ITERS::Int64 = 20, POW_TOL = 1e-8) where{T}
    (alpha <= 0 || alpha >= 1) && throw(DomainError("The exponent α of the power cone has to be in (0, 1)."))
    new(3, alpha, MAX_ITERS, POW_TOL)
  end
end
PowerCone(args...) = PowerCone{DefaultFloat}(args...)


# The projection onto the power cone is described in
# Hien - Differential properties of Euclidean projections onto power cone (2015)
function project!(v::AbstractVector{T}, cone::PowerCone{T}) where{T}
  # Check the special cases first
  # 1. v in K_pow => v = v
  in_cone(v, cone, 0.) && return nothing

  # 2. -v in K_pow^* => v .= 0
  if in_dual(-v, cone, 0.)
    v .= zero(T)
    return nothing
  end

  # 3. v not in K_pow and -v not in K_pow^* and z == 0 => x = max(x, 0), y = max(y, 0)
  if abs(v[3]) <= cone.POW_TOL
    v[1] = max(v[1], 0)
    v[2] = max(v[2], 0)
    return nothing
  end

  # 4. Otherwise solve the following problem
  # find  r
  # s.t.  sigma(x0, y0, z0, r) == 0
  # and   0 < r < ||z0||
  #
  # x = 0.5 * (x0 + sqrt(x0^2 + 4α*r*(||z0||-r)))
  # y = 0.5 * (y0 + sqrt(y0^2 + 4(1-α)*r*(||z0||-r)))
  # z = z0 * r / ||z0||
  project_pow!(v, cone)
end

 # find the zero of above condition for r by applying Newton's method
  function project_pow!(v::AbstractVector{T}, cone::PowerCone{T}) where{T}
  x0 = v[1]
  y0 = v[2]
  z0 = v[3]
  r = abs(z0) / 2
  ϕx = zero(T)
  ϕy = zero(T)
  # compute a zero of phi(v, r, α)
  for k = 1:cone.MAX_ITER
    ϕx = ϕc(x0, z0, r, cone.α)
    ϕy = ϕc(y0, z0, r, 1 - cone.α)
    phi = ϕ(ϕx, ϕy, r, cone.α)
    abs(phi) < cone.POW_TOL && break

    dϕx_dr = dϕc_dr(ϕx, x0, z0, r, cone.α)
    dϕy_dr = dϕc_dr(ϕy, y0, z0, r, 1 - cone.α)
    dphi_dr = dϕ_dr(ϕx, ϕy, dϕx_dr, dϕy_dr, r, cone.α)
    r = r - phi / dphi_dr

    # ensure 0 < r < abs(z0)
    r = min(max(r, 0), abs(z0))
  end

  # given a solution r, update the vector components (x, y, z)
  v[1] = ϕx
  v[2] = ϕy
  v[3] = z0 * r / abs(z0)
  return nothing
end

function ϕc(x0::T, z0::T, r::T, α::T) where{T <: Real}
  return max(0.5 * (x0 + sqrt(x0^2 + 4 * α * r * (abs(z0) - r))), 1e-10)
end

function dϕc_dr(ϕx::T, x0::T, z0::T, r::T, α::T) where{T <: Real}
  return α / (2 * ϕx - x0) * (abs(z0) - 2 * r)
end

function ϕ(ϕx::T, ϕy::T, r::T, α::T) where{T <: Real}
  return ϕx^α * ϕy^(1 - α) - r
end

# dϕ / dr = Π fi^αi * (Σ αi fi' / fi) - 1
function dϕ_dr(ϕx::T, ϕy::T, ϕx_dr::T, ϕy_dr::T, r::T, α::T) where{T <: Real}
  return ϕx^α * ϕy^(1-α) * (α * ϕx_dr / ϕx + (1 - α) * ϕy_dr / ϕy) - 1
end

function in_cone(v::AbstractVector{T}, cone::PowerCone{T}, tol::T) where{T}
  x = v[1]
  y = v[2]
  z = v[3]
  α = cone.α
  return x >= 0 && y >= 0 && x^α * y^(1 - α) >= abs(z) - tol
end

# Kpow^* = { (s, t, w) | (s / α)^α * (t/(1-α))^(1-α) >= abs(w), s >= 0, t >= 0 }
function in_dual(v::AbstractVector{T}, cone::PowerCone{T}, tol::T) where{T}
  s = v[1]
  t = v[2]
  w = v[3]
  α = cone.α
  return s >= -tol && t >= -tol && s^α * t^(1-α) >= abs(w) * α^α * (1-α)^(1-α) - tol
end

function in_pol_recc(v::AbstractVector{T},cone::PowerCone{T}, tol::T) where{T}
  return in_dual(-v, cone, tol)
end


"""
    DualExponentialCone(MAX_ITERS::Int64 = 100, EXP_TOL = 1e-8)

Creates the dual exponential cone ``\\mathcal{K}^*_{exp} = \\{(x, y, z) \\mid x < 0,  -xe^{y/x} \\leq e^1 z \\} \\cup \\{ (0,y,z) \\mid   y \\geq 0, z \\geq 0 \\}``
"""
struct DualExponentialCone{T} <: AbstractConvexCone{T}
  dim::Int
  v0::Vector{T}
  primal_cone::ExponentialCone{T}

  function DualExponentialCone{T}(dim::Int64 = 3, MAX_ITERS::Int64 = 100, EXP_TOL = 1e-8) where{T}
    new(3, zeros(T, 3), ExponentialCone{T}(dim, MAX_ITERS, EXP_TOL))
  end
end
DualExponentialCone(args...) = DualExponentialCone{DefaultFloat}(args...)

"""
    DualPowerCone(alpha::Float64, MAX_ITERS::Int64 = 20, POW_TOL = 1e-8)

Creates the 3-d dual power cone ``\\mathcal{K}^*_{pow} = \\{(u, v, w) \\mid \\left( \\frac{u}{\\alpha}\\right)^\\alpha \\left( \\frac{v}{1-\\alpha}\\right)^{(1-\\alpha)} \\geq  \\|w\\|, u \\geq 0, v \\geq 0 \\}`` with ``0 < \\alpha < 1``
"""
struct DualPowerCone{T} <: AbstractConvexCone{T}
  dim::Int
  v0::Vector{T}
  primal_cone::PowerCone{T}

  function DualPowerCone{T}(alpha::Float64, MAX_ITERS::Int64 = 20, POW_TOL = 1e-8) where{T}
    (alpha <= 0 || alpha >= 1) && throw(DomainError, "The exponent α of the dual power cone has to be in (0, 1).")
    new(3, zeros(T,3), PowerCone{T}(alpha, MAX_ITERS, POW_TOL))
  end
end
DualPowerCone(args...) = DualPowerCone{DefaultFloat}(args...)
DualCones = Union{DualExponentialCone, DualPowerCone}
in_cone(v::AbstractVector{T}, cone::DualCones, tol::Real) where {T <: Real} = in_dual(v, cone.primal_cone, tol)
in_dual(v::AbstractVector{T}, cone::DualCones, tol::Real) where {T <: Real} = in_cone(v, cone.primal_cone, tol)
in_pol_recc(v::AbstractVector{T}, cone::DualCones, tol::Real) where {T <: Real} = in_dual(-v, cone, tol)

# Project dual cones by using Moreau decomposition: Proj^*(v) = v + Proj(-v)
function project!(v::AbstractVector{T}, cone::DualCones) where{T}
  @. cone.v0 = v
  @. v *= -1.0
  project!(v, cone.primal_cone)
  @. v += cone.v0
end

# Union for all types where the user has to provide extra information to create the cone
ArgumentCones = Union{PowerCone, DualPowerCone}


# ----------------------------------------------------
# Box
# ----------------------------------------------------
"""
    Box(l, u)

Creates a box or intervall with lower boundary vector ``l \\in  \\mathbb{R}^m \\cup \\{-\\infty\\}^m`` and upper boundary vector``u \\in \\mathbb{R}^m\\cup \\{+\\infty\\}^m``.
"""
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
support_function!(x::AbstractVector{T}, B::Box{T}, tol::T) where{T} = support_function(x, B, tol)

function in_pol_recc(x::AbstractVector{T}, B::Box{T}, tol::T) where{T}
    !any(XU -> (XU[2] == Inf && XU[1] > tol), zip(x,B.u)) && !any(XL -> (XL[2] == -Inf && XL[1] < -tol), zip(x,B.l))
end

function scale!(box::Box{T}, e::AbstractVector{T}) where{T}
	@. box.l = box.l * e
	@. box.u = box.u * e
	return nothing
end

function Base.deepcopy(box::Box{T}) where {T}
  Box{T}(deepcopy(box.l), deepcopy(box.u))
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

function support_function!(x::SplitVector{T}, C::CompositeConvexSet{T}, tol::T) where{T}
	sum(xC -> support_function!(xC[1], xC[2], tol), zip(x.views, C.sets))
end

function in_pol_recc(x::SplitVector{T}, C::CompositeConvexSet{T}, tol::T) where{T}
	all(xC -> in_pol_recc(xC[1], xC[2], tol), zip(x.views, C.sets))
end

function in_pol_recc!(x::SplitVector{T}, C::CompositeConvexSet{T}, tol::T) where{T}
  all(xC -> in_pol_recc!(xC[1], xC[2], tol), zip(x.views, C.sets))
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

# An in-place method that is faster, but uses the input variable y as workspace
function support_function!(y::SplitView{T}, cone::AbstractConvexCone{T}, tol::T) where{T}
  @. y *= - one(T)
  in_dual!(y, cone, tol) ? 0. : Inf;
end

# Notice: for every convex set apart from PsdCone and PsdConeTriangle use the normal non-modifying function
# for PsdCone and PsdConeTriangle we are using (faster) in-place functions.
function in_dual!(x::AbstractVector{T}, cone::AbstractConvexSet{T}, tol::T) where{T}
  return in_dual(x, cone, tol)
end

function in_pol_recc!(x::AbstractVector{T}, cone::AbstractConvexSet{T}, tol::T) where{T}
  return in_pol_recc(x, cone, tol)
end

function scale!(cone::AbstractConvexCone{T}, ::AbstractVector{T}) where{T}
  return nothing
end

function rectify_scaling!(E, work, set::Union{SecondOrderCone{T}, PsdCone{T}, DensePsdCone{T}, PsdConeTriangle{T}, PsdConeTriangleLOBPCG{T}, DensePsdConeTriangle{T}, PowerCone{T}, DualPowerCone{T}, ExponentialCone{T}, DualExponentialCone{T}}) where{T}
  return rectify_scalar_scaling!(E, work)
end
rectify_scaling!(E, work, set::Union{ZeroSet{<:Real}, Nonnegatives{<:Real}, Box{<:Real}}) = false

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
