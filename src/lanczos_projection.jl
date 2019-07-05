using Arpack, LinearMaps, LinearAlgebra
# ----------------------------------------------------
# Positive Semidefinite Cone
# ----------------------------------------------------
mutable struct PsdConeTriangleLanczos{T} <: AbstractConvexCone{T}
  dim::Int
  n::Int
  positive_subspace::Bool
  X::Matrix{T}  # Matrix under projection
  Z::Matrix{T}  # Ritz vectors
  λ::Vector{T}  # Ritz values
  λ_rem::T
  z_rem::Vector{T}
  buffer_size::Int
  iter_number::Int
  # History
  residual_history::Vector{T}
  λ_rem_history::Vector{T}
  subspace_dim_history::Vector{Int}
  λ_rem_multiplications::Vector{Int}

  function PsdConeTriangleLanczos{T}(dim::Int) where{T}
      dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
      n = Int(1/2*(sqrt(8*dim + 1) - 1)) # Solution of (n^2 + n)/2 = length(x) obtained by WolframAlpha
      n*(n + 1)/2 == dim || throw(DomainError(dim, "dimension must be a square"))
      new(dim, n, true,
        zeros(T, n, n), # X
        zeros(T, n, 0), # Ζ
        zeros(T, n), # λ
        zero(T), # λ_rem
        randn(T, n), # z_rem
        3, # buffer_size
        0, # iter_number
        zeros(T, 0), # residual_history
        zeros(T, 0), # λ_rem_history
        zeros(Int, 0), # subspace_dim_history
        zeros(Int, 0) # λ_rem_multiplications
        )
  end
end
PsdConeTriangleLanczos(dim) = PsdConeTriangleLanczos{DefaultFloat}(dim)

function project_to_nullspace!(x::AbstractVector{T}, tmp::AbstractVector{T}, U::Array{T}) where {T}
  # Project x to the nullspace of U', i.e. x .= (I - U*U')*x
  # tmp is a vector used in intermmediate calculations
  mul!(tmp, U', x)
  BLAS.gemv!('N', -one(T), U, tmp, one(T), x)
end

function estimate_λ_rem(X::Symmetric{T, Matrix{T}}, U::Matrix{T}, n::Int, x0::Vector{T}) where T
  # Estimates largest eigenvalue of the Symmetric X on the subspace we discarded
  # Careful, we need to enforce all iterates to be orthogonal to the range of U
  tmp = zeros(T, size(U, 2))
  offset = T(100)
  function custom_mul!(y::AbstractVector{T}, x::AbstractVector{T}) where {T}
      # Performs y .= (I - U*U')*X*x
      # y .= X*x - U*(U'*(X*x))
      project_to_nullspace!(x, tmp, U)
      mul!(y, X, x)
      project_to_nullspace!(y, tmp, U)
      axpy!(offset, x, y)
  end
  project_to_nullspace!(x0, tmp, U)
  A = LinearMap{T}(custom_mul!, size(X, 1); ismutating=true, issymmetric=true)
  λ_rem, v_rem, nconv, niter, nmult, resid = eigs(A, nev=n, ncv=20, which=:LR, tol=1e-4, v0=x0)
  return λ_rem .- offset, v_rem, nmult
end

function generate_subspace(X::Symmetric{T, Matrix{T}}, cone::PsdConeTriangleLanczos{T}) where {T}
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

function generate_subspace_chol(X::Symmetric{T, Matrix{T}}, cone::PsdConeTriangleLanczos{T}) where {T}
  XZ = X*cone.Z; ZXZ = cone.Z'*XZ
  V, XV, VXV, VX2V = propagate_subspace(cone.Z, XZ, ZXZ, cone.λ)
  if size(VXV, 2) == 0
      return cone.Z, XZ
  end
  R = cholesky([I VXV; VXV' VX2V], Val(true)) #check=false)
  #=
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
  =#
  @show logdet(R)
  W = (R.L\[V XV]')';
  @show norm(W'*W - I)
  XW = X*W

  return W, XW
end

function project!(x::AbstractArray, cone::PsdConeTriangleLanczos{T}) where{T}
  n = cone.n
  cone.iter_number += 1

  if mod(cone.iter_number, 40) == 0 || size(cone.Z, 2) == 0 || size(cone.Z, 2) >= cone.n/2 || n == 1
      append!(cone.λ_rem_multiplications, 0)
      return project_exact!(x, cone)
  end

  if !cone.positive_subspace
      @. x = -x
  end
  populate_upper_triangle!(cone.X, x)
  X = Symmetric(cone.X)

  # W, XW = generate_subspace(X, cone)
  W, XW = generate_subspace_chol(X, cone)
  Xsmall = W'*XW
  l, V, first_positive, first_negative = eigen_sorted(Symmetric(Xsmall), 1e-10);

  # Positive Ritz pairs
  Vp = V[:, first_positive:end]
  U = W*Vp; λ = l[first_positive:end];

  buffer_idx = max(first_negative-cone.buffer_size-1,1):max(first_negative,0)
  Ub = W*V[:, buffer_idx]; λb = l[buffer_idx]

  # Projection
  Xπ = Symmetric(U*Diagonal(λ)*U');
  
  # Residual Calculation
  # R = XW*Vp - U*Diagonal(λ)
  R = X*U - U*Diagonal(λ)

  # Important: why not [U Ub]?
  λ_rem, z_rem, nmult = estimate_λ_rem(X, U, 1, cone.z_rem)
  append!(cone.λ_rem_multiplications, nmult)
  append!(cone.λ_rem_history, maximum(λ_rem))
  λ_rem .= max.(λ_rem, 0.0)
  eig_sum = sum(λ_rem).^2 + (n - size(W, 2) - length(λ_rem))*minimum(λ_rem).^2
  residual = sqrt(2*norm(R)^2 + eig_sum)
  #=
  @show norm(R, 2)
  @show sqrt(eig_sum)
  =#
  
  if !cone.positive_subspace
    @inbounds for j in 1:n, i in 1:j
      Xπ.data[i, j] -= X[i, j]
    end
  end
  extract_upper_triangle!(Xπ, x)

  # It should be [U ub z_rem]
  cone.Z = [U z_rem]
  # cone.z_rem = randn(size(U, 1))
  cone.λ = [λ; λ_rem]

  append!(cone.residual_history, norm(R))
  append!(cone.subspace_dim_history, size(cone.Z, 2))
end

function project_exact!(x::AbstractArray{T}, cone::PsdConeTriangleLanczos{T}) where{T}
  n = cone.n

  # handle 1D case
  if length(x) == 1
      x = max.(x,zero(T))
  else
      # symmetrized square view of x
      X = cone.X
      populate_upper_triangle!(X, x)

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
      extract_upper_triangle!(X, x)
      
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
          idx = findfirst(λ[sorted_idx] .< 0) # First negative index
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
  append!(cone.residual_history, 0.0)
  append!(cone.λ_rem_history, 0.0)
  append!(cone.subspace_dim_history, size(cone.Z, 2))
  return nothing
end

function eigen_sorted(A::Symmetric, tol::AbstractFloat=0.0)
  λ, U = eigen(A)
  sorted_idx = sortperm(λ)
  λ = λ[sorted_idx]
  U = U[:, sorted_idx]
  first_positive = findfirst(λ[sorted_idx] .> tol)
  if isa(first_positive, Nothing)
      first_positive = length(λ) + 1
  end
  first_negative = findfirst(λ[sorted_idx] .< tol)
  if isa(first_negative, Nothing)
      first_negative = 0
  end

  return λ, U, first_positive, first_negative
end