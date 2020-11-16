export NoRegularizer, TikonovRegularizer, FrobeniusNormRegularizer, Type1, Type2, RollingMemory, RestartedMemory, AndersonAccelerator, EmptyAccelerator

# An abstract type for fixed-point acceleration methods
# Fixed point problem x = g(x) with residual f(x) = x - g(x)
"""
    AbstractAccelerator{T} where {T <: AbstractFloat}

Abstract supertype for acceleration objects that can be used to speed up a fixed-point iterations g = g(x) of a nonexpansive operator `g`. They must implement the following functions:
  - update_history!(aa::AbstractAccelerator{T}, g::AbstractVector{T}, x::AbstractVector{T})
  - accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AbstractAccelerator, num_iter)
"""
abstract type AbstractAccelerator end

get_mem(:: AbstractAccelerator) = 0
is_active(::AbstractAccelerator) = true
# ---------------------------
# AndersonAccelerator
# ---------------------------

"""
    AbstractRegularizer

Abstract supertype for Anderson Acceleration Type-II regularisation schemes.
"""
abstract type AbstractRegularizer end
struct TikonovRegularizer <: AbstractRegularizer end
struct FrobeniusNormRegularizer <: AbstractRegularizer end
struct NoRegularizer <: AbstractRegularizer end

"""
    AbstractBroydenType

Abstract supertype for the Broyden type of the accelerator.
"""
abstract type AbstractBroydenType end
struct Type1 <: AbstractBroydenType end
struct Type2 <: AbstractBroydenType end

"""
    AbstractMemory

Abstract supertype for the memory management of the accelerator.
"""
abstract type AbstractMemory end

"""
    RollingMemory

The accelerator will append new iterates to the history and discard the oldest iterate.
"""
struct RollingMemory <: AbstractMemory end

"""
    RestartedMemory

The accelerator will delete the history once the memory buffers are full.
"""
struct RestartedMemory <: AbstractMemory end


"""
    AndersonAccelerator{T, R, BT, M} <: AbstractAccelerator{T}

Accelerator object implementing Anderson Acceleration. Parameterized by:

 - T: AbstractFloat, floating-point type
 - R: AbstractRegularizer
 - BT: Broyden-type, i.e. Type-I or Type-II
 - M: AbstractMemory, how full memory buffers are handled
"""
mutable struct AndersonAccelerator{T, R, BT, M} <: AbstractAccelerator
  init_phase::Bool
  mem::Int64
  dim::Int64
  iter::Int64
  num_accelerated_steps::Int64
  fail_counter::Array{Int64}
  fail_eta::Array{Int64}
  fail_singular::Array{Int64}
  cond::Array{Float64, 2}
  x_last::AbstractVector{T}
  g_last::AbstractVector{T}
  f::AbstractVector{T}
  f_last::AbstractVector{T}
  eta::AbstractVector{T}
  F::AbstractMatrix{T}
  X::AbstractMatrix{T}
  G::AbstractMatrix{T}
  M::AbstractMatrix{T}
  λ::T # regularisation parameter
  start_iter::Int64 # start acceleration after x iterations of the main algorithm
  start_accuracy::T # start acceleration after required accuracy
  activated::Bool # a flag that shows that the accelerator has been started
  reg_log::Vector{T} # log regularisation size
  update_time::Float64
  accelerate_time::Float64

  function AndersonAccelerator{T, R, BT, M}() where {T <: AbstractFloat, R <: AbstractRegularizer, BT <: AbstractBroydenType, M <: AbstractMemory}
    new(true, 0, 0, 0, 0, zeros(Int64, 0), zeros(Int64, 0), zeros(Int64, 0), zeros(T, 0, 2), zeros(T, 1), zeros(T, 1),  zeros(T, 1), zeros(T, 1), zeros(T, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zero(T), 0, 0., false, zeros(T, 0), 0., 0.)
  end

  function AndersonAccelerator{T, R, BT, M}(dim::Int64; mem::Int64 = 5, λ = 1e-8, start_iter::Int64 = 2, start_accuracy::T = Inf) where {T <: AbstractFloat, R <: AbstractRegularizer, BT <: AbstractBroydenType, M <: AbstractMemory}
    mem <= 2 && throw(DomainError(mem, "Memory has to be bigger than two."))
    dim <= 0 && throw(DomainError(dim, "Dimension has to be a positive integer."))

    # mem shouldn't be bigger than the dimension
    mem = min(mem, dim)
    new(true, mem, dim, 0, 0, zeros(Int64,0), zeros(Int64,0), zeros(Int64,0), zeros(Float64, 0, 2), zeros(T,dim), zeros(T, dim), zeros(T, dim),  zeros(T, dim), zeros(T, mem), zeros(T, dim, mem), zeros(T, dim, mem), zeros(T, dim, mem), zeros(T, mem, mem), λ, start_iter, start_accuracy, false, zeros(T, 0), 0., 0.)
  end
end
# define some default constructors for parameters
AndersonAccelerator(args...; kwargs...)  = AndersonAccelerator{Float64, NoRegularizer, Type2, RollingMemory}(args...; kwargs...)
AndersonAccelerator{T}(args...; kwargs...) where {T <: AbstractFloat} = AndersonAccelerator{T, NoRegularizer, Type2, RollingMemory}(args...; kwargs...)
AndersonAccelerator{BT}(args...; kwargs...) where {BT <: AbstractBroydenType} = AndersonAccelerator{Float64, NoRegularizer, BT, RollingMemory}(args...; kwargs...)
AndersonAccelerator{M}(args...; kwargs...) where {M <: AbstractMemory} = AndersonAccelerator{Float64, NoRegularizer, Type2, M}(args...; kwargs...)
AndersonAccelerator{R}(args...; kwargs...) where {R <: AbstractRegularizer} = AndersonAccelerator{Float64, R, Type2, RollingMemory}(args...; kwargs...)

get_type(::AndersonAccelerator{T, R, BT, M}) where {T,R, BT, M} = BT
get_memory(::AndersonAccelerator{T, R, BT, M}) where {T,R, BT, M} = M
get_regularizer(::AndersonAccelerator{T, R, BT, M}) where {T,R, BT, M} = R
get_mem(aa::AndersonAccelerator) = aa.mem 


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
  # accelerator_time = 0. ?
end


function empty_caches!(aa::AndersonAccelerator)
  aa.F .= 0;
  aa.X .= 0;
  aa.G .= 0;
  aa.iter = 0 #this is important as for the RestartedMemory it holds information how many rows are full
end

is_active(aa::AndersonAccelerator) = aa.active

function check_activation!(aa::AndersonAccelerator, num_iter::Int64)
  if aa.active
    return nothing
  else
    if num_iter >= aa.start_iter
      aa.active = true
    end
  end
  return nothing
end
check_activation!(aa::AbstractAccelerator, num_iter::Int64) = nothing
check_activation!(aa::AbstractAccelerator, r_prim::T, r_dual::T) where {T <: AbstractFloat} = nothing

function check_activation!(aa::AndersonAccelerator, r_prim::T, r_dual::T) where {T <: AbstractFloat}
  if aa.active
    return nothing
  else
    if r_prim <= aa.start_accuracy && r_dual <= aa.start_iter
      aa.active = true
    end
  end
  return nothing
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


"""
  update_history!(aa, g, x)

- Update history of accelerator `acc` with iterates g = g(xi) and xi
- Computes residuals f = x - g
"""
function update_history!(aa::AndersonAccelerator{T, R, BT, M}, g::AbstractVector{T}, x::AbstractVector{T}) where {T <: AbstractFloat, R <: AbstractRegularizer, BT <: AbstractBroydenType, M <: AbstractMemory}
  if aa.init_phase
    @. aa.x_last = x
    @. aa.g_last = g
    @. aa.f_last = x - g
    aa.init_phase = false
    return nothing
  end

  update_time_start = time()


  j = (aa.iter % aa.mem) + 1 # (aa.iter % aa.mem) number of cols filled, j is the next col where data should be entered

  if j == 1 && aa.iter != 0
    apply_memory_approach!(aa) # for a RestartedMemory approach we want to flush the data cache matrices and start from scratch
  end

  # compute residual
  @. aa.f = x - g

  # fill memory with deltas
  @. aa.X[:, j] = x - aa.x_last
  @. aa.G[:, j] = g - aa.g_last
  @. aa.F[:, j] = aa.f - aa.f_last

  # depending on method type
  assemble_inv_matrix!(aa, aa.X, aa.F)

  # set previous values for next iteration
  @. aa.x_last = x
  @. aa.g_last = g
  @. aa.f_last = aa.f

  aa.iter += 1
  aa.update_time += time() - update_time_start
  return nothing
end

apply_memory_approach!(aa::AndersonAccelerator{T, R, BT, RollingMemory}) where {T, R, BT} = false
function apply_memory_approach!(aa::AndersonAccelerator{T, R, BT, RestartedMemory}) where {T, R, BT}
    empty_caches!(aa)
    return true
end

function assemble_inv_matrix!(aa::AndersonAccelerator{T, R, Type1}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat, R <: AbstractRegularizer}
  aa.M[:, :] = aa.X' * aa.F
end

function assemble_inv_matrix!(aa::AndersonAccelerator{T, R, Type2}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat, R <: AbstractRegularizer}
  aa.M[:, :] = aa.F' * aa.F
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


function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T, R}, num_iter) where {T <: AbstractFloat, R <: AbstractRegularizer}

  l = min(aa.iter, aa.mem) #number of columns filled with data
  l < 3 && return true
  accelerate_time_start = time()

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

  initialise_eta!(eta, aa, X, F)
  # aa.cond = vcat(aa.cond, [num_iter cond(M)])

  # aa.eta = aa.M  \ (X' * f) (type1)
  info = solve_linear_sys!(M, X, F, eta, aa)

  if (info < 0 || norm(eta, 2) > 1e4)
    if info < 0
      push!(aa.fail_singular, num_iter)
      #@warn("Acceleration failed at aa.iter: $(aa.iter) because of info.")
    elseif norm(eta, 2) > 1e4
      push!(aa.fail_eta, num_iter)
      # @warn("Acceleration failed at aa.iter: $(aa.iter) because of norm(eta).")
    end
    push!(aa.fail_counter, num_iter)
    aa.accelerate_time += time() - accelerate_time_start
    return false
  else
     # @show(eta)
    aa.num_accelerated_steps += 1
    g[:] = g - G * eta
    aa.accelerate_time += time() - accelerate_time_start
    return true
  end
end

function initialise_eta!(eta::AbstractVector{T}, aa::AndersonAccelerator{T, R, Type1}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat, R <: AbstractRegularizer}
  eta[:] = X' * aa.f
end

function initialise_eta!(eta::AbstractVector{T}, aa::AndersonAccelerator{T, R, Type2}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat, R <: AbstractRegularizer}
  eta[:] = F' * aa.f
end

function solve_linear_sys!(M::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, R}) where {T <: AbstractFloat, R <: TikonovRegularizer}
  # add regularisation term
  for i = 1:size(M, 1)
    M[i, i] += aa.λ
  end
  # solve regularised problem
  info = _gesv!(M, eta)
end

function solve_linear_sys!(M::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, R}) where {T <: AbstractFloat, R <: FrobeniusNormRegularizer}
  β = aa.λ * (norm(X)^2 + norm(F)^2)
  push!(aa.reg_log, β)
  # add regularisation term
  for i = 1:size(M, 1)
    M[i, i] += β
  end
  # solve regularised problem
  info = _gesv!(M, eta)
end

function solve_linear_sys!(M::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, R}) where {T <: AbstractFloat, R <: NoRegularizer}
  info = _gesv!(M, eta)
end



function print_failure_rate(aa::AndersonAccelerator{<: AbstractFloat})
  println("AA - Failure rate: Failed = $(aa.fail_counter), Total iter = $(aa.iter - 1) ($(round(aa.fail_counter / (aa.iter - 1) * 100, digits = 2)) %)")
end


# ---------------------------
# EmptyAccelerator
# ---------------------------

struct EmptyAccelerator{T} <: AbstractAccelerator 
  function EmptyAccelerator{T}() where {T <: AbstractFloat}
      return new{T}()
  end
end
EmptyAccelerator(args...; kwargs...) = EmptyAccelerator{Float64}(args...; kwargs...)
EmptyAccelerator{T}(dim::Int64) where {T <: AbstractFloat} = EmptyAccelerator{T}()

function update_history!(ea::EmptyAccelerator{T}, x::AbstractVector{T}, g::AbstractVector{T}) where {T <: AbstractFloat}
  return nothing
end

function empty_history!(ea::EmptyAccelerator{T}) where {T <: AbstractFloat}
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::EmptyAccelerator{T}, iter) where {T <: AbstractFloat}
  return true
end
