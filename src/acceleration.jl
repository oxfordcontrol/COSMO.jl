export NoRegularizer, TikonovRegularizer, FrobeniusNormRegularizer, Type1, Type2, RollingMemory, RestartedMemory, AndersonAccelerator, EmptyAccelerator
export IterActivation, AccuracyActivation, IterOrAccuracyActivation, ImmediateActivation
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
is_actived(::AbstractAccelerator) = true
is_safeguarding(::AbstractAccelerator) = false

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

abstract type AbstractActivationReason end
"Activate accelerator after `start_iter` iterations of the main algorithm."
struct IterActivation <: AbstractActivationReason
  start_iter::Int64
  function IterActivation(start_iter::Int64)
    start_iter < 2 && ArgumentError("start_iter has to be at least 2.")
    new(start_iter)
  end
end
"Activate accelerator after accuracy of main algorithm <= `start_accuracy`."
struct AccuracyActivation <: AbstractActivationReason
  start_accuracy::Float64
  function AccuracyActivation(start_accuracy::Float64)
    start_accuracy >= 0. && ArgumentError("start_accuracy has to be a non-negative number.")
    new(start_accuracy)
  end 
end
"Activate accelerator after accuracy of main algorithm <= `start_accuracy` or iter >= `start_iter`."
struct IterOrAccuracyActivation <: AbstractActivationReason
  start_accuracy::Float64
  start_iter::Int64
  function IterOrAccuracyActivation(start_accuracy::Float64, start_iter::Int64)
    start_accuracy >= 0. && ArgumentError("start_accuracy has to be a non-negative number.")
    start_iter < 2 && ArgumentError("start_iter has to be at least 2.")
    new(start_accuracy, start_iter)
  end 
end
"Activate accelerator immediately."
struct ImmediateActivation <: AbstractActivationReason end



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
  w_acc::AbstractVector{T}
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
  τ::T # safeguarding slack parameters  
  activation_reason::AbstractActivationReason 
  activated::Bool # a flag that shows that the accelerator has been started
  safeguarded::Bool 
  reg_log::Vector{T} # log regularisation size
  update_time::Float64
  accelerate_time::Float64

  function AndersonAccelerator{T, R, BT, M}() where {T <: AbstractFloat, R <: AbstractRegularizer, BT <: AbstractBroydenType, M <: AbstractMemory}
    new(true, 0, 0, 0, 0, zeros(Int64, 0), zeros(Int64, 0), zeros(Int64, 0), zeros(T, 0, 2), zeros(T, 1), zeros(T, 1), zeros(T, 1),  zeros(T, 1), zeros(T, 1), zeros(T, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zero(T), T(1.1), ImmediateActivation(), false, false, zeros(T, 0), 0., 0.)
  end

  function AndersonAccelerator{T, R, BT, M}(dim::Int64; mem::Int64 = 5, λ = 1e-8, start_iter::Int64 = 2, start_accuracy::T = Inf, safeguarded::Bool = false, τ::T = 1.1, activation_reason::AbstractActivationReason = ImmediateActivation()) where {T <: AbstractFloat, R <: AbstractRegularizer, BT <: AbstractBroydenType, M <: AbstractMemory}
    mem <= 2 && throw(DomainError(mem, "Memory has to be bigger than two."))
    dim <= 0 && throw(DomainError(dim, "Dimension has to be a positive integer."))

    # mem shouldn't be bigger than the dimension
    mem = min(mem, dim)
    new(true, mem, dim, 0, 0, zeros(Int64,0), zeros(Int64,0), zeros(Int64,0), zeros(Float64, 0, 2), zeros(T, dim), zeros(T,dim), zeros(T, dim), zeros(T, dim),  zeros(T, dim), zeros(T, mem), zeros(T, dim, mem), zeros(T, dim, mem), zeros(T, dim, mem), zeros(T, mem, mem), λ, τ, activation_reason, false, safeguarded, zeros(T, 0), 0., 0.)
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
is_safeguarding(aa::AndersonAccelerator) = aa.safeguarded

function empty_history!(aa::AndersonAccelerator{T}) where {T <: AbstractFloat}
  aa.F .= 0;
  aa.X .= 0;
  aa.G .= 0;

  # aa.f .= 0; we need it for safeguarding
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

is_actived(aa::AndersonAccelerator) = aa.activated


# dispatch on the activation reason
check_activation!(aa::AndersonAccelerator, args...; kwargs...) = check_activation!(aa, aa.activation_reason, args...; kwargs...)

function check_activation!(aa::AndersonAccelerator, activation_reason::ImmediateActivation, num_iter::Int64, r_prim::T, r_dual::T) where {T <: AbstractFloat}
    (!aa.activated && num_iter >= 2) && (aa.activated = true)
end

function check_activation!(aa::AndersonAccelerator, activation_reason::IterActivation, num_iter::Int64, r_prim::T, r_dual::T) where {T <: AbstractFloat}
  (!aa.activated && num_iter >= activation_reason.start_iter) && (aa.activated = true)
end

function check_activation!(aa::AndersonAccelerator, activation_reason::AccuracyActivation, num_iter::Int64, r_prim::T, r_dual::T) where {T <: AbstractFloat}
  (!aa.activated && r_prim <= activation_reason.start_accuracy && r_dual <= activation_reason.start_accuracy) && (aa.activated = true)
end

function check_activation!(aa::AndersonAccelerator, activation_reason::IterOrAccuracyActivation, num_iter::Int64, r_prim::T, r_dual::T) where {T <: AbstractFloat}
  if !aa.activated
    if r_prim <= activation_reason.start_accuracy && r_dual <= activation_reason.start_accuracy
      aa.activated = true
    elseif num_iter >= activation_reason.start_iter
      aa.activated = true
    end
  end
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


function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T, R}, num_iter; rws::Union{Nothing, ResidualWorkspace} = nothing, ws::Union{Workspace{T}, Nothing} = nothing   ) where {T <: AbstractFloat, R <: AbstractRegularizer}

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
    elseif norm(eta, 2) > 1e4
      push!(aa.fail_eta, num_iter)
    end
    push!(aa.fail_counter, num_iter)
    aa.accelerate_time += time() - accelerate_time_start
    return false
  else
    aa.num_accelerated_steps += 1
    # calculate the accelerated candidate point
    @. aa.w_acc = g 
    aa.w_acc -= G * eta
    aa.accelerate_time += time() - accelerate_time_start
  	# safeguard the acceleration
    if aa.safeguarded
      nrm_f_acc = fixed_point_residual_norm(rws, ws, aa.w_acc)
      @show(nrm_f_acc, norm(aa.f, 2))
        if nrm_f_acc <= aa.τ * norm(aa.f, 2)  #acc.f = (w_prev - w)
          @. g = aa.w_acc
          println("Point accepted")
        else
          println("Point declined")
        end
    else # or just overwrite anyway
      @. g = aa.w_acc	
    end
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

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::EmptyAccelerator{T}, args...; kwargs... ) where {T <: AbstractFloat}
  return true
end
check_activation!(aa::EmptyAccelerator, args...; kwargs...) = nothing
