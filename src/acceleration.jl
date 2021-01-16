export NoRegularizer, TikonovRegularizer, FrobeniusNormRegularizer, Type1, Type2, RollingMemory, RestartedMemory, AndersonAccelerator, EmptyAccelerator
export QRDecomp, NormalEquations
# An abstract type for fixed-point acceleration methods
# Fixed point problem x = g(x) with residual f(x) = x - g(x)
"""
    AbstractAccelerator

Abstract supertype for acceleration objects that can be used to speed up a fixed-point iterations g = g(x) of a nonexpansive operator `g`. 
They must implement the following functions to communicate with the fixed-point algorithm:
  - update_history!(aa::AbstractAccelerator, g, x, num_iter) #stores the fixed-point iterates
  - accelerate!(g::AbstractVector, x, aa::AbstractAccelerator, num_iter) #recombines past iterates to determine an accelerated point and overwrites `g`
  - restart!(aa::AbstractAccelerator, args...; kwargs...) # algorithm wants the accelerator to restart
  - activate!(aa::AbstractAccelerator, args...; kwargs...) # algorithm activates the accelerator
  - log!(aa::AbstractAccelerator, args...; kwargs...) # algorithm tells accelerator to log certain information for debugging

The algorithm has to be able to query the following information:
  - was_successful(aa::AbstractAccelerator) #indicate whether accelerate! was succesful at the last iteration
  - get_memory_size(aa::AbstractAccelerator) #return the memory length
  - is_active(aa::AbstractAccelerator) #returns whether the accelerator is active
"""
abstract type AbstractAccelerator end

# update_history!(aa::AbstractAccelerator, g, x, iter) = nothing
# accelerate!(g, x, aa::AbstractAccelerator, iter) = nothing
# restart!(aa::AbstractAccelerator, args...; kwargs...) = nothing
# get_memory_size(aa::AbstractAccelerator) = 0
# is_active(aa::AbstractAccelerator) = true
# was_successful(aa::AbstractAccelerator) = false
# activate!(aa::AbstractAccelerator) = nothing
log!(aa::AbstractAccelerator, args...; kwargs...) = false

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
    AbstractLeastSquaresMethod

Abstract supertype for least squares solution method when Type-II acceleration is used
"""
abstract type AbstractLeastSquaresMethod end
struct NormalEquations <: AbstractLeastSquaresMethod end
struct QRDecomp <: AbstractLeastSquaresMethod end

"""
    AbstractBroydenType

Abstract supertype for the Broyden type of the accelerator.
"""
abstract type AbstractBroydenType end
struct Type1 <: AbstractBroydenType end
struct Type2{LSM} <: AbstractBroydenType 
  function Type2{LSM}() where {LSM <: AbstractLeastSquaresMethod}
    return new{LSM}()
  end
end

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
    AndersonAccelerator{T, RE, BT, MT} <: AbstractAccelerator

Accelerator object implementing Anderson Acceleration. Parameterized by:

 - T: AbstractFloat, floating-point type
 - RE: AbstractRegularizer
 - BT: Broyden-type, i.e. Type1 or Type2
 - MT: AbstractMemory, how full memory buffers are handled
"""
mutable struct AndersonAccelerator{T, RE, BT, MT} <: AbstractAccelerator
  init_phase::Bool
  mem::Int64
  min_mem::Int64
  dim::Int64
  iter::Int64
  num_accelerated_steps::Int64
  fail_counter::Array{Int64}
  fail_eta::Array{Int64}
  fail_singular::Array{Int64}
  w_acc::Vector{T}
  x_last::Vector{T}
  g_last::Vector{T}
  f::Vector{T}
  f_last::Vector{T}
  eta::Vector{T}
  F::Matrix{T}
  X::Matrix{T}
  G::Matrix{T}
  M::Matrix{T}
  Q::Matrix{T}
  R::Matrix{T}
  λ::T # regularisation parameter
  activated::Bool # a flag that shows that the accelerator has been started
  success::Bool # a flag to indicate whether the last attempted acceleration was successful
  reg_log::Vector{T} # log regularisation size
  # logging 
  update_time::Float64
  accelerate_time::Float64
  acc_post_time::Float64
  acceleration_status::Vector{Tuple{Int64, Symbol}}
  safeguarding_status::Vector{Tuple{Int64, T, T, T}} 
  num_safe_accepted::Int64
  num_safe_declined::Int64
  activate_logging::Bool
  
  function AndersonAccelerator{T, RE, BT, MT}() where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
    new(true, 0, 3, 0, 0, 0, zeros(Int64, 0), zeros(Int64, 0), zeros(Int64, 0), zeros(T, 1), zeros(T, 1), zeros(T, 1),  zeros(T, 1), zeros(T, 1), zeros(T, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zeros(T, 1, 1), zero(T), false, false, zeros(T, 0), 0., 0., 0., Vector{Tuple{Int64, Symbol}}(undef, 0),  Vector{Tuple{Int64, T, T, T}}(undef, 0),  0, 0, false)
  end

  function AndersonAccelerator{T, RE, BT, MT}(dim::Int64; mem::Int64 = 10, min_mem::Int64 = 3, λ::T = T(1e-8)) where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
    mem <= 1 && throw(DomainError(mem, "Memory has to be bigger than one."))
    dim <= 0 && throw(DomainError(dim, "Dimension has to be a positive integer."))

    # mem shouldn't be bigger than the dimension
    mem = min(mem, dim)
    min_mem = min(min_mem, mem)

    F = zeros(T, dim, mem)
    X = zeros(T, dim, mem)
    G = zeros(T, dim, mem)
    M = zeros(T, mem, mem)
    Q = zeros(T, 0, 0)
    R = zeros(T, 0, 0)
    x_last = zeros(T, dim) 
    new(true, mem, min_mem, dim, 0, 0, zeros(Int64,0), zeros(Int64,0), zeros(Int64,0), x_last, zeros(T,dim), zeros(T, dim), zeros(T, dim),  zeros(T, dim), zeros(T, mem), F, X, G, M, Q, R, λ, false, false, zeros(T, 0), 0., 0., 0., Vector{Tuple{Int64, Symbol}}(undef, 0), Vector{Tuple{Int64, T, T, T}}(undef, 0), 0, 0, false)
  end

  function AndersonAccelerator{T, RE, Type2{QRDecomp}, MT}(dim::Int64; mem::Int64 = 10, min_mem::Int64 = 3, λ::T = T(1e-8)) where {T <: AbstractFloat, RE <: AbstractRegularizer, MT <: AbstractMemory}
    mem <= 1 && throw(DomainError(mem, "Memory has to be bigger than one."))
    dim <= 0 && throw(DomainError(dim, "Dimension has to be a positive integer."))

    # mem shouldn't be bigger than the dimension
    mem = min(mem, dim)
    min_mem = min(min_mem, mem)

    # for QR decomposition we don't need some of the caches
    x_last = zeros(T, 0)
    F = zeros(T, 0, 0)
    X = zeros(T, 0, 0)
    G = zeros(T, dim, mem)
    M = zeros(T, 0, 0)
    Q = zeros(T, dim, mem)
    R = zeros(T, mem, mem)
  
    new(true, mem, min_mem, dim, 0, 0, zeros(Int64,0), zeros(Int64,0), zeros(Int64,0), x_last, zeros(T,dim), zeros(T, dim), zeros(T, dim),  zeros(T, dim), zeros(T, mem), F, X, G, M, Q, R, λ, false, false, zeros(T, 0), 0., 0., 0., Vector{Tuple{Int64, Symbol}}(undef, 0), Vector{Tuple{Int64, T, T, T}}(undef, 0), 0, 0, false)
  end

end
# define some default constructors for parameters
AndersonAccelerator(args...; kwargs...)  = AndersonAccelerator{Float64, NoRegularizer, Type2{QRDecomp}, RestartedMemory}(args...; kwargs...)
AndersonAccelerator{T}(args...; kwargs...) where {T <: AbstractFloat} = AndersonAccelerator{T, NoRegularizer, Type2{QRDecomp}, RestartedMemory}(args...; kwargs...)
AndersonAccelerator{BT}(args...; kwargs...) where {BT <: AbstractBroydenType} = AndersonAccelerator{Float64, NoRegularizer, BT, RollingMemory}(args...; kwargs...)
AndersonAccelerator{M}(args...; kwargs...) where {M <: AbstractMemory} = AndersonAccelerator{Float64, NoRegularizer, Type2{NormalEquations}, M}(args...; kwargs...)
AndersonAccelerator{R}(args...; kwargs...) where {R <: AbstractRegularizer} = AndersonAccelerator{Float64, R, Type2{NormalEquations}, RollingMemory}(args...; kwargs...)


get_type(::AndersonAccelerator{T, R, BT, M}) where {T,R, BT, M} = BT
get_memory(::AndersonAccelerator{T, R, BT, M}) where {T,R, BT, M} = M
get_regularizer(::AndersonAccelerator{T, R, BT, M}) where {T,R, BT, M} = R
get_memory_size(aa::AndersonAccelerator) = aa.mem 
is_active(aa::AndersonAccelerator) = aa.activated
was_successful(aa::AndersonAccelerator) = aa.success


function log!(aa::AndersonAccelerator, iter::Int64, status::Symbol)
  aa.activate_logging && push!(aa.acceleration_status, (num_iter, status))
  if status == :acc_guarded_declined
    aa.num_safe_declined += 1
  elseif status == :acc_guarded_accepted
    aa.num_accelerated_steps += 1
    aa.num_safe_accepted += 1
  elseif status == :acc_unguarded
    aa.num_accelerated_steps += 1   
  end
  return nothing
end


function activate!(aa::AndersonAccelerator) 
  aa.activated = true
end

function restart!(aa::AndersonAccelerator, iter, status)
  empty_history!(aa)
  log!(aa, iter, status)
end



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
 
end

"Reset the pointer that tracks the number of valid cols in the cached history of past vectors."
function empty_caches!(aa::AndersonAccelerator)
  #to save time we can leave the old data in here, this is a potential source of problems if the indexing of valid data gets messed up
  # aa.F .= 0; 
  # aa.X .= 0;
  # aa.G .= 0;
  aa.iter = 0 #this is important as for the RestartedMemory holds information on how many rows are full
end

"""
  update_history!(aa, g, x)

- Update history of accelerator `aa` with iterates g = g(xi)
- Computes residuals f = x - g
"""
function update_history!(aa::AndersonAccelerator{T, RE, BT, MT}, g::AbstractVector{T}, x::AbstractVector{T}, num_iter::Int64) where {T <: AbstractFloat, RE <: AbstractRegularizer, BT <: AbstractBroydenType, MT <: AbstractMemory}
  update_time_start = time()  
  if aa.activated
    # compute residual
    @. aa.f = x - g

    if aa.init_phase
      set_prev_vectors!(aa.x_last, aa.g_last, aa.f_last, x, g, aa.f)
      aa.init_phase = false
      return nothing
    end

    j = (aa.iter % aa.mem) + 1 # (aa.iter % aa.mem) number of cols filled, j is the next col where data should be entered

    if j == 1 && aa.iter != 0
      apply_memory_approach!(aa, num_iter) # for a RestartedMemory approach we want to flush the data cache matrices and start from scratch
    end

    # store Δx, Δg, Δf in X, G, F
    fill_caches!(j, aa.X, aa.G, aa.F, x, g, aa.f, aa.x_last, aa.g_last, aa.f_last)

    # set previous values for next iteration
    set_prev_vectors!(aa.x_last, aa.g_last, aa.f_last, x, g, aa.f)

    aa.iter += 1
    aa.update_time += time() - update_time_start
  end
  return nothing
end

# This method is a copy of the method above to test TypeII with QR decomposition
function update_history!(aa::AndersonAccelerator{T, RE, Type2{QRDecomp}, RestartedMemory}, g::AbstractVector{T}, x::AbstractVector{T}, num_iter::Int64) where {T <: AbstractFloat, RE <: AbstractRegularizer}
  if aa.activated
    update_time_start = time()  

    # compute residual
    @. aa.f = x - g

    if aa.init_phase
      set_prev_vectors!(aa.g_last, aa.f_last, g, aa.f)
      aa.init_phase = false
      return nothing
    end

    j = (aa.iter % aa.mem) + 1 # (aa.iter % aa.mem) number of cols filled, j is the next col where data should be entered

    if j == 1 && aa.iter != 0
      apply_memory_approach!(aa, num_iter) # for a RestartedMemory approach we want to flush the data cache matrices and start from scratch
    end

    #G[:, j] = Δg
    fill_delta!(j, aa.G, g, aa.g_last)
    
    # use f_last memory to store Δf = f - f_last 
    compute_Δf!(aa.f_last, aa.f)

    # QR decomposition step
    qr!(aa.Q, aa.R, aa.f_last, j)

    # set previous values for next iteration
    set_prev_vectors!(aa.g_last, aa.f_last, g, aa.f)

    aa.iter += 1
    aa.update_time += time() - update_time_start
  end
  return nothing
end

# adapted from https://github.com/JuliaNLSolvers/NLsolve.jl
# add new column to QR-Factorisation: F = QR 
function qr!(Q::AbstractMatrix{T}, R::AbstractMatrix{T}, Δf::AbstractVector{T}, j::Int64) where {T <: AbstractFloat}
  
  n, m = size(Q)
  n == length(Δf) || throw(DimensionMismatch())
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ j ≤ m || throw(ArgumentError())

  @inbounds for k in 1:(j - 1)
    qk = uview(Q, :, k)
    ip = dot(qk, Δf)
    R[k, j] = ip

    axpy!(-ip, qk, Δf)
  end
  @inbounds begin
    nrm_f = norm(Δf, 2)
    R[j, j] = nrm_f
    @. Q[:, j] = Δf / nrm_f
  end
  return nothing
end


" Update the history matrices X = [Δxi, Δxi+1, ...], G = [Δgi, Δgi+1, ...] and F = [Δfi, Δfi+1, ...]."
function fill_caches!(j::Int64, X::AbstractMatrix{T}, G::AbstractMatrix{T}, F::AbstractMatrix{T}, x::AbstractVector{T}, g::AbstractVector{T}, f::AbstractVector{T}, x_last::AbstractVector{T}, g_last::AbstractVector{T}, f_last::AbstractVector{T}) where {T <: AbstractFloat}
  # fill memory with deltas
  @. X[:, j] = x - x_last # Δx
  @. G[:, j] = g - g_last # Δg
  @. F[:, j] = f - f_last # Δf
  return nothing
end

" Update a single history matrices V = [vgi, vgi+1, ...] ."
function fill_delta!(j::Int64, V::AbstractMatrix{T}, v::AbstractVector{T}, v_last::AbstractVector{T}) where {T <: AbstractFloat}
  @inbounds for i in eachindex(v)
    V[i, j] = v[i] - v_last[i]
  end
end

"Compute Δf = f - f_last and store result in f_last."
function compute_Δf!(f_last::AbstractVector{T}, f::AbstractVector{T}) where {T <: AbstractFloat}
  
  @inbounds for i in eachindex(f_last)
    f_last[i] *= -one(T)
    f_last[i] += f[i]
  end
end

"Store a copy of the last x, g, f to be able to compute Δx, Δg, Δf at the next step."
function set_prev_vectors!(x_last::AbstractVector{T}, g_last::AbstractVector{T}, f_last::AbstractVector{T}, x::AbstractVector{T}, g::AbstractVector{T}, f::AbstractVector{T}) where {T <: AbstractFloat}
  @. x_last = x
  @. g_last = g
  @. f_last = f
  return nothing
end

function set_prev_vectors!(g_last::AbstractVector{T}, f_last::AbstractVector{T}, g::AbstractVector{T}, f::AbstractVector{T}) where {T <: AbstractFloat}
  @. g_last = g
  @. f_last = f
  return nothing
end

"Depending on the AndersonAccelerator parameter AbstractMemory, dispatch on the correct method that handles the case when the memory buffers are full."
function apply_memory_approach!(aa::AndersonAccelerator{T, R, BT, RestartedMemory}, num_iter) where {T, R, BT}
    empty_caches!(aa)
    aa.activate_logging && log_restart!(aa, num_iter, :memory_full)
    return true
end
apply_memory_approach!(aa::AndersonAccelerator{T, R, BT, RollingMemory}, num_iter) where {T, R, BT} = nothing

function assemble_inv_matrix!(W::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, aa::AndersonAccelerator{T, R, Type1, M}) where {T <: AbstractFloat, R <: AbstractRegularizer, M <: AbstractMemory}
  mul!(W, X', F)
end

function assemble_inv_matrix!(W::AbstractMatrix{T}, X::AbstractMatrix{T}, F::AbstractMatrix{T}, aa::AndersonAccelerator{T, R, Type2{NormalEquations}, M}) where {T <: AbstractFloat, R <: AbstractRegularizer, M <: AbstractMemory}
  mul!(W, F', F)
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


function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T, RE}, num_iter::Int64) where {T <: AbstractFloat, RE <: AbstractRegularizer}
  aa.success = false
  l = min(aa.iter, aa.mem) #number of columns filled with data
  if l < aa.min_mem
     aa.activate_logging && push!(aa.acceleration_status, (num_iter, :not_enough_cols))
     return nothing
  end
  accelerate_time_start = time()

  if l < aa.mem
    eta = uview(aa.eta, 1:l)
    X = uview(aa.X, :, 1:l)
    M = uview(aa.M, 1:l, 1:l)
    G = uview(aa.G, :, 1:l)
    F = uview(aa.F, :, 1:l)
  else
    eta = aa.eta
    X = aa.X
    M = aa.M
    G = aa.G
    F = aa.F
  end
  assemble_inv_matrix!(M, X, F, aa)

  # depending on method type
  initialise_eta!(eta, aa, X, F)

  # aa.eta = aa.M  \ (X' * f) (type1)
  info = solve_linear_sys!(M, X, F, eta, aa)

  # num_iter < 100 && num_iter > 0 && @show(num_iter, eta)
  # if (info < 0 || norm(eta, 2) > 1e4)
    # if info < 0
      # push!(aa.fail_singular, num_iter)
      # aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_singular))
    if norm(eta, 2) > 1e4
      # push!(aa.fail_eta, num_iter)
      aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_eta_norm))
    # end
    # push!(aa.fail_counter, num_iter)
    
    aa.accelerate_time += time() - accelerate_time_start
    return nothing
  else
    # calculate the accelerated candidate point
    # g = g - G * eta
    BLAS.gemv!('N', -one(T), G, eta, one(T), g)
    aa.accelerate_time += time() - accelerate_time_start
    aa.success = true
    return nothing
  end
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::AndersonAccelerator{T, RE, Type2{QRDecomp}}, num_iter::Int64) where {T <: AbstractFloat, RE <: AbstractRegularizer}
  aa.success = false
  l = min(aa.iter, aa.mem) #number of columns filled with data
  if l < aa.min_mem
     aa.activate_logging && push!(aa.acceleration_status, (num_iter, :not_enough_cols))
     return nothing
  end
  accelerate_time_start = time()

  eta = uview(aa.eta, 1:l)
  G = uview(aa.G, :, 1:l)
  Q = uview(aa.Q, :, 1:l)
  R = uview(aa.R, 1:l, 1:l) 

  # solve least squares problem ||f_k - η Fk ||_2 where Fk = QR 
  # initialise_eta!(eta, aa, X, F)
  mul!(eta, Q', aa.f) # vec(aa.f)? 
  info = solve_linear_sys!(R, eta, aa)
  if info < 0
    aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_cond_r))
    aa.accelerate_time += time() - accelerate_time_start
    return nothing
  end
  # num_iter < 100 && num_iter > 0 && @show(num_iter, eta)
  
  # TODO: maybe replace this with a check of the condition number of R
  if norm(eta, 2) > 1e4
    aa.activate_logging && push!(aa.acceleration_status, (num_iter, :fail_eta_norm))
    aa.accelerate_time += time() - accelerate_time_start
    return nothing
  else
    # calculate the accelerated candidate point
    # g = g - G * eta
    BLAS.gemv!('N', -one(T), G, eta, one(T), g)
        
    aa.accelerate_time += time() - accelerate_time_start
    aa.success = true
    return nothing
  end
end



function initialise_eta!(eta::AbstractVector{T}, aa::AndersonAccelerator{T, R, Type1}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat, R <: AbstractRegularizer}
  mul!(eta, X', aa.f)
end

function initialise_eta!(eta::AbstractVector{T}, aa::AndersonAccelerator{T, R, Type2{NormalEquations}}, X::AbstractMatrix{T}, F::AbstractMatrix{T}) where {T <: AbstractFloat, R <: AbstractRegularizer}
  mul!(eta, F', aa.f)
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


function solve_linear_sys!(R::AbstractMatrix{T}, eta::AbstractVector{T}, aa::AndersonAccelerator{T, Re, Type2{QRDecomp}}) where {T <: AbstractFloat, Re <: NoRegularizer}
  try
    LAPACK.trtrs!('U', 'N', 'N', R, eta)
    # LinearAlgebra.ldiv!(R, eta)
    # LinearAlgebra.BLAS.trsv!('U', 'N', 'N', R, eta) # seems to be equally fast
    return 1
   catch
     return -1
  end
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

function update_history!(ea::EmptyAccelerator{T}, args...; kwargs...) where {T <: AbstractFloat}
  return nothing
end

function empty_history!(ea::EmptyAccelerator{T}) where {T <: AbstractFloat}
  return nothing
end

function accelerate!(g::AbstractVector{T}, x::AbstractVector{T}, aa::EmptyAccelerator{T}, args...; kwargs... ) where {T <: AbstractFloat}
  return true
end
was_successful(aa::EmptyAccelerator) = false
log_restart!(aa::AbstractAccelerator, iter::Int64, reason::Symbol) = nothing