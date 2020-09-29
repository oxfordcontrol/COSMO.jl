# An abstract type for fixed-point acceleration methods
# Fixed point problem x = g(x) with residual f(x) = x - g(x)
abstract type AbstractAccelerator{T <: AbstractFloat}  end

mutable struct AndersonAccelerator{T} <: AbstractAccelerator{T}
  is_type1::Bool
  init_phase::Bool
  mem::Int64
  dim::Int64
  iter::Int64
  fail_counter::Array{Int64}
  cond::Float64
  x_last::AbstractVector{T}
  g_last::AbstractVector{T}
  f::AbstractVector{T}
  f_last::AbstractVector{T}
  eta::AbstractVector{T}
  F::AbstractMatrix{T}
  X::AbstractMatrix{T}
  G::AbstractMatrix{T}
  M::AbstractMatrix{T}

  function AndersonAccelerator{T}() where {T <: AbstractFloat}
    new(true, true, 0, 0, 0, zeros(Int64,0), 0., zeros(T, 1), zeros(T, 1), zeros(T, 1),  zeros(T, 1), zeros(T, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1))
  end

  function AndersonAccelerator{T}(dim::Int64; mem::Int64 = 4, is_type1::Bool = true) where {T <: AbstractFloat}
    mem <= 0 && throw(DomainError(mem, "Memory has to be a positive integer."))
    dim <= 0 && throw(DomainError(dim, "Dimension has to be a positive integer."))
    new(is_type1, true, mem, dim, 0, zeros(Int64,0), 0., zeros(T,dim), zeros(T, dim), zeros(T, dim),  zeros(T, dim), zeros(T, mem), zeros(T, dim, mem), zeros(T, dim, mem), zeros(T, dim, mem), zeros(T, mem, mem))
  end

end

function empty_history!(aa::AndersonAccelerator{T}) where {T <: AbstractFloat}
  aa.F .= 0;
  aa.X .= 0;
  aa.G .= 0;

  aa.f .= 0;
  aa.f_last .= 0;
  aa.g_last .= 0;
  aa.x_last .= 0;
  aa.eta .= 0;

  aa.iter = 0
  aa.init_phase = true
  aa.cond = 0.
end


function compute_alphas(eta::Vector{T}) where {T <: AbstractFloat}
  n = length(eta)
  alpha = zeros(T, n+1)
  alpha[1] = eta[1]
  for i = 2:n
    alpha[i] = eta[i] - eta[i-1]
  end
  alpha[n + 1] = 1 - eta[n]
  return alpha
end


function update_history!(aa::AbstractAccelerator{T}, g::AbstractVector{T}, x::AbstractVector{T}) where {T <: AbstractFloat}
  if aa.init_phase
    @. aa.x_last = x
    @. aa.g_last = g
    @. aa.f_last = x - g
    aa.init_phase = false
    return nothing
  end
  j = (aa.iter % aa.mem) + 1

  # compute residual
  @. aa.f = x - g

  # fill memory with deltas
  @. aa.X[:, j] = x - aa.x_last
  @. aa.G[:, j] = g - aa.g_last
  @. aa.F[:, j] = aa.f - aa.f_last

  if aa.is_type1
    aa.M[:, :] = aa.X' * aa.F
  else
    aa.M[:, :] = aa.F' * aa.F
  end
  aa.cond = cond(aa.M)

  # set previous values for next iteration
  @. aa.x_last = x
  @. aa.g_last = g
  @. aa.f_last = aa.f

  aa.iter += 1
end

# BLAS gesv! wrapper with error handling
# solve A X = B and for X and save result in B
function _gesv!(A, B)
  try
    LinearAlgebra.LAPACK.gesv!(A, B)
    return 1
   catch
     return -1
  end
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T}, num_iter) where {T <: AbstractFloat}
  l = min(aa.iter, aa.mem)
  l < 2 && return true

  # aa.iter > 100 && return true
  if l < aa.mem
    eta = view(aa.eta, 1:l)
    X = view(aa.X, :, 1:l)
    M = view(aa.M, 1:l, 1:l)
    G = view(aa.G, :, 1:l)
    F = view(aa.F, :, 1:l)
  else
    eta = aa.eta
    X = aa.X
    M = aa.M
    G = aa.G
    F = aa.F
  end

  if aa.is_type1
    eta[:] = X' * aa.f
  else
    eta[:] = F' * aa.f
  end
  # aa.eta = aa.M  \ (X' * f) (type1)
  info = _gesv!(M, eta)

  if (info < 0 || norm(eta, 2) > 1e4)
    #@warn("Acceleration failed at aa.iter: $(aa.iter)")
    push!(aa.fail_counter, num_iter)
    return false
  else
    g[:] = g - G * eta
    return true
  end
end

function print_failure_rate(aa::AndersonAccelerator{<: AbstractFloat})
  println("AA - Failure rate: Failed = $(aa.fail_counter), Total iter = $(aa.iter - 1) ($(round(aa.fail_counter / (aa.iter - 1) * 100, digits = 2)) %)")
end

struct EmptyAccelerator{T} <: AbstractAccelerator{T} end

function update_history!(ea::EmptyAccelerator{T}, x::AbstractVector{T}, g::AbstractVector{T}) where {T <: AbstractFloat}
  return nothing
end

function empty_history!(ea::EmptyAccelerator{T}) where {T <: AbstractFloat}
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::EmptyAccelerator{T}, iter) where {T <: AbstractFloat}
  return true
end
