using Arpack, LinearMaps, LinearAlgebra

# ----------------------------------------------------
# Positive Semidefinite Cone
# ----------------------------------------------------
mutable struct PsdConeTriangleLanczos{T} <: AbstractConvexCone{T}
  dim::Int
  n::Int
  iter_number::Int
  positive_subspace::Bool
  X::Matrix{T}
  data::LOBPCGIterable{T}
  # History
  residual_history::Vector{T}
  λ_rem_history::Vector{T}
  subspace_dim_history::Vector{Int}
  λ_rem_multiplications::Vector{Int}

  function PsdConeTriangleLanczos{T}(dim::Int) where{T}
      dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
      n = Int(1/2*(sqrt(8*dim + 1) - 1)) # Solution of (n^2 + n)/2 = length(x) obtained by WolframAlpha
      n*(n + 1)/2 == dim || throw(DomainError(dim, "dimension must be a square"))
      initial_dim = min(10, Int(floor(n/2)))
      data = LOBPCGIterable(zeros(T, n, n))
      data.X = randn(T, n, initial_dim)
      new(dim, n, 0, true,
        zeros(T, n, n), # X
        data,
        zeros(T, 0), # residual_history
        zeros(T, 0), # λ_rem_history
        zeros(Int, 0), # subspace_dim_history
        zeros(Int, 0) # λ_rem_multiplications
        )
  end
end
PsdConeTriangleLanczos(dim) = PsdConeTriangleLanczos{DefaultFloat}(dim)

function get_tolerance(cone::PsdConeTriangleLanczos{T}) where T
  return T(sqrt(cone.n)/cone.iter_number^(1.05))*10
end

function project!(x::AbstractArray, cone::PsdConeTriangleLanczos{T}) where{T}
  n = cone.n
  cone.iter_number += 1

  if !cone.positive_subspace
      @. x = -x
  end
  populate_upper_triangle!(cone.X, x)
  tol = get_tolerance(cone)
  # @show cone.X
  try
    initialize!(cone.data, -cone.X, cone.data.X)
    cone.data.tol = tol
    cone.data, status = lobpcg!(cone.data, 10)
    cone.data.λ = -cone.data.λ
    if status != :converged
      if !cone.positive_subspace
        # revert flipping
        @. x = -x
      end
      append!(cone.λ_rem_multiplications, 0)
      return project_exact!(x, cone)
    end
    # lobpcg!(cone.data, tol)
  catch e
    if isa(e, ArgumentError) || isa(e, PosDefException)
      # @warn "Exact projection"
      # Perform "exact" projection
      if !cone.positive_subspace
        # revert flipping
        @. x = -x
      end
      append!(cone.λ_rem_multiplications, 0)
      return project_exact!(x, cone)
    else
      throw(e)
    end 
  end
  # return project_exact!(x, cone)

  append!(cone.residual_history, 0)
  append!(cone.subspace_dim_history, size(cone.data.X, 2))

  
  # Reconstruct projection
  U = copy(cone.data.X)
  rmul!(U, Diagonal(sqrt.(max.(cone.data.λ, 0))))
  if cone.positive_subspace
    BLAS.syrk!('U', 'N', 1.0, U, 0.0, cone.X)
  else
    BLAS.syrk!('U', 'N', 1.0, U, -1.0, cone.X)
  end
  keep = min(sum(cone.data.λ .> 1e-6) + max(Int(floor(n/50)), 3), length(cone.data.λ))
  cone.data.X = cone.data.X[:, 1:keep]

  extract_upper_triangle!(cone.X, x)
end

function set_memory()

end

function project_exact!(x::AbstractArray{T}, cone::PsdConeTriangleLanczos{T}) where{T}
  n = cone.n

  # handle 1D case
  if length(x) == 1
      x = max.(x,zero(T))
  else
      # symmetrized square view of x
      X = zeros(T, n, n)
      populate_upper_triangle!(cone.X, x)

      # compute eigenvalue decomposition
      # then round eigs up and rebuild
      λ, U  = eigen!(Symmetric(cone.X))
      Up = U[:, λ .> 0]
      sqrt_λp = sqrt.(λ[λ .> 0])
      rmul!(Up, Diagonal(sqrt_λp))
      BLAS.syrk!('U', 'N', 1.0, Up, 0.0, cone.X)
      extract_upper_triangle!(cone.X, x)
      
      # Save the subspace that we will be tracking
      if sum(λ .> 0) <= sum(λ .< 0)
        cone.positive_subspace = true
      else
        λ .= -λ # Equivalent to considering -cone.X instead of cone.X
        cone.positive_subspace = false
      end
      sorted_idx = sortperm(λ)
      idx = findfirst(λ[sorted_idx] .> 0) # First positive index
      if isa(idx, Nothing)
          idx = length(λ) + 1
      end
      # Take also a few vectors from the discarted eigenspace
      idx = max(idx - 5, 1)
      cone.data.X = U[:, sorted_idx[idx:end]]
      cone.data.m = size(cone.data.X, 2)
      cone.data.λ = λ[sorted_idx[idx:end]]
  end
  append!(cone.residual_history, 0.0)
  append!(cone.λ_rem_history, 0.0)
  append!(cone.subspace_dim_history, size(cone.data.X, 2))
  #@show "Done with Exact projection", size(cone.Z.Q1, 2), cone.n
  # cone.z_rem .= randn(n)
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