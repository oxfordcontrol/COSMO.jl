using LinearAlgebra, Statistics
using Printf
using Random

mutable struct LOBPCGIterable{T}
    n::Int
    m::Int
    m_max::Int
    iteration::Int
    A::Matrix{T}
    X::Matrix{T}
    R::Matrix{T}
    P::Matrix{T}

    AX::Matrix{T}
    AP::Matrix{T}

    λ::Vector{T}
    condition_estimates::Vector{T}
    active_indices::BitVector
    tol::T
    verbosity::Int
    is_restarted::Bool
    buffer_size::Int
    largest::Bool

    _memoryGramm::Matrix{T}
    _memoryR::Matrix{T}
    _memoryP::Matrix{T}
    _memoryAR::Matrix{T}
    _memoryAP::Matrix{T}

    function LOBPCGIterable(A::Matrix{T};
        tol::T=T(1e-4), verbosity=1,
        buffer_size::Int=max(Int(floor(size(A, 1) / 50)), 3),
        max_dim=Int(ceil(size(A, 1)/4)),
        largest=true
        ) where T
        n = size(A, 1)
        @assert max_dim <= Int(ceil(n/4)) "Too many (max) requested eigenvalues. Use exact eigendecomposition or use exact eigendecomposition"

        new{T}(n, 0, max_dim, 0, A,
            Matrix{T}(undef, n, 0), # X
            Matrix{T}(undef, n, max_dim), Matrix{T}(undef, n, 0), # R, P
            Matrix{T}(undef, n, 0), Matrix{T}(undef, n, 0), # ΑΧ, AP
            zeros(T, 0), # λ
            [-log10(eps(T))/2], # condition_estimates
            BitArray(undef, 0), # active_indices
            tol, verbosity,
            true, # is_restarted
            buffer_size, # buffer
            largest,
            Matrix{T}(undef, 3*max_dim, 6*max_dim), # _memoryGramm
            Matrix{T}(undef, n, max_dim), Matrix{T}(undef, n, max_dim), # _memoryR/P
            Matrix{T}(undef, n, max_dim), Matrix{T}(undef, n, max_dim)  # _memoryAR/AP
        )
    end
end

function initialize!(data::LOBPCGIterable{Float64}, A::Matrix{Float64}, X::Matrix{Float64})
    @assert data.n == size(A, 1) == size(A, 2)
    data.A = A
    initialize!(data, X)
end

function initialize!(data::LOBPCGIterable{T}, X::Matrix{T}) where T
    data.m = size(X, 2)
    @assert size(X, 1) == data.n
    data.m > data.m_max && throw(ArgumentError("Too many columns in initial matrix - try increasing max_dim."))

    data.X = X
    data.A = data.A
    data.AX = Symmetric(data.A)*X
    data.λ, _ = RayleighRitz!(data.X, data.AX)
    data.active_indices = BitArray(undef, data.m)
    data.active_indices .= true

    data.P = similar(X)
    data.AP = similar(X)

    data.iteration = 0
    data.condition_estimates = [-log10(eps(T))/2]
    restart!(data)
end

function lobpcg!(data::LOBPCGIterable{Float64}, max_iter=10)
    flag = :not_converged
    while data.iteration < max_iter && flag == :not_converged
        data.iteration += 1
        flag = iterate!(data)
        if data.verbosity > 2
            @printf("Iteration %i current block size %i \n",
                        data.iteration, sum(data.active_indices));
            println("Eigenvalues:", data.λ)
            # println("Residuals  :", data.residual_norms)
        end
    end
    if flag != :converged
        data.verbosity > 1 && @warn string(sum(data.active_indices), " eigenpairs not converged!")
    end

    return data, flag
end

function lobpcg(A::Matrix{Float64}, X::Matrix{Float64}, max_iter=10; kwargs...)
    data = LOBPCGIterable(A; kwargs...)
    initialize!(data, A, X)
    lobpcg!(data, max_iter)
end

function iterate!(data::LOBPCGIterable{T}) where T
    # data.R = data.AX - data.X*Diagonal(data.λ)
    # copyto!(data.R, data.AX)
    data.R = copy(data.AX)
    @inbounds @simd for j = 1:data.m
        axpy!(-data.λ[j], view(data.X, :, j), view(data.R, :, j))
    end
    residual_norms = column_norms(data.R)

    R, P, AR, AP, residual_norms = drop_converged_indices!(data, residual_norms)
    # WARNING: AR is uninitialized
    if length(R) == 0
        return :converged
    end
    # R -= data.X*(data.X'*R)
    BLAS.gemm!('N', 'N', -one(T), data.X, data.X' * R, one(T), R)
    
    _, success = orthonormalize!(R)
    # @assert size(R, 2) - size(AR, 2) == 0
    if !success
        data.verbosity > 1 && @warn "Ill conditioned residual. Recomputing."
        R = Matrix(qr(R).Q)
        explicit_computation = true
    else
        explicit_computation = false
    end
    ### Main computational time in the following time
    mul!(AR, Symmetric(data.A), R) # Symmetric might not offer much when data.m >> 1
    ###

    if !data.is_restarted
        F, success = orthonormalize!(P)
        if success
            rdiv!(AP, F.U)
        else
            data.verbosity > 1 && @warn "Orthonormalization failed. Restarting"
            P, AP = restart!(data)
        end
    end
    
    condition_estimate = zero(T)
    explicit_computation |= all(residual_norms .<= 3 * eps(T) / 5)
    G, Q = compute_gramm_matrix!(data, explicit_computation, R, P, AR, AP)
    for condition_try = 1:2
        condition_estimate = log10(cond(Q)) + 1;
        condition_estimate_mean = mean(data.condition_estimates[
            max(1, data.iteration - 10 - Int(round(log(size(R, 2))))):data.iteration])
    
        if (condition_estimate / condition_estimate_mean > 2 && condition_estimate > 2) || condition_estimate > 8
            data.verbosity > 1 && @info string("Restart on iteration: ", data.iteration, " due to:") condition_estimate
            if !data.is_restarted && condition_try == 1
                P, AP = restart!(data)
                G, Q = compute_gramm_matrix!(data, explicit_computation, R, P, AR, AP)
            else
                data.verbosity > 0 && @warn "Gramm matrix ill-conditioned: results unpredictable"
            end
        else
            break
        end
    end
    push!(data.condition_estimates, condition_estimate)

    λ, W = eigen!(G, Q) # 2nd most computationally intensive line
    flag = update_subspaces!(data, λ, W, R, AR, P, AP)
    flag && return :excessive_columns
    
    return :not_converged
end

function restart!(data::LOBPCGIterable{T}) where T
    data.is_restarted = true
    P = view(data._memoryP, :, 1:0)
    AP = view(data._memoryAP, :, 1:0)
    return P, AP
end

function update_subspaces!(data::LOBPCGIterable{T}, λ, W, R, AR, P, AP) where T
    if data.largest
        m_new = min(max(sum(λ .> 0) + data.buffer_size, data.m), length(λ))
        data.λ = λ[end-m_new+1:end]
        W = view(W, :, size(W, 2) - m_new + 1:size(W,2))
    else
        m_new = min(max(sum(λ .< 0) + data.buffer_size, data.m), length(λ))
        data.λ = λ[1:m_new]
        W = view(W, :, 1:m_new)
    end
    m_new > min(data.n/4, data.m_max) && return true

    data.P = [R P] * W[data.m+1:end, :]
    data.X = data.X * W[1:data.m, :] + data.P
    data.AP = [AR AP] * W[data.m+1:end, :]
    data.AX = data.AX * W[1:data.m, :] + data.AP
    data.is_restarted = false
    if data.verbosity > 1 && m_new > data.m
        @info string("Increasing number of eigenvalues to ", m_new)
        restart!(data)
    end
    data.m = m_new
    return false
end

function compute_gramm_matrix!(data::LOBPCGIterable{T}, explicit_gram::Bool, R, P, AR, AP) where T
    # Computes and returns the following two (symmetric) matrices (memory) efficiently
    # G = Symmetric([data.X data.R data.P]' * data.A * [data.X data.R data.P])
    # Q = Symmetric([data.X data.R data.P]' * [data.X data.R data.P])

    gram_dimension = data.m + size(R, 2) + size(P, 2)
    G = view(data._memoryGramm, 1:gram_dimension, 1:gram_dimension)
    G11, G12, G13, G22, G23, G33 = get_gramm_views(data, G, R, P)

    if !explicit_gram
        set_to_diagonal!(G11, data.λ)
    else
        mul!(G11, data.X', data.AX)
    end
    mul!(G12, data.X', AR)
    mul!(G13, data.X', AP)
    mul!(G22, R', AR)
    mul!(G23, R', AP)
    mul!(G33, P', AP)

    Q = view(data._memoryGramm, 1:gram_dimension, 3*data.m+1:3*data.m+gram_dimension)
    Q11, Q12, Q13, Q22, Q23, Q33 = get_gramm_views(data, Q, R, P)

    if !explicit_gram
        set_to_diagonal!(Q11, one(T))
        set_to_diagonal!(Q22, one(T))
        Q12 .= 0
    else
        BLAS.syrk!('U', 'T', one(T), data.X, zero(T), Q11)
        BLAS.syrk!('U', 'T', one(T), R, zero(T), Q22)
        mul!(Q12, data.X', R)
    end
    BLAS.syrk!('U', 'T', one(T), P, zero(T), Q33)
    mul!(Q13, data.X', P)
    mul!(Q23, R', P)

    return Symmetric(G), Symmetric(Q)
end

function get_gramm_views(data, A, R, P)
    m, k, l = data.m, size(R, 2), size(P, 2)
    A11 = view(A, 1:m, 1:m)
    A12 = view(A, 1:m, m+1:m+k)
    A13 = view(A, 1:m, m+k+1:m+k+l)
    A22 = view(A, m+1:m+k, m+1:m+k)
    A23 = view(A, m+1:m+k, m+k+1:m+k+l)
    A33 = view(A, m+k+1:m+k+l, m+k+1:m+k+l)
    return A11, A12, A13, A22, A23, A33
end

function set_to_diagonal!(A::AbstractArray, d::Vector{T}) where {T}
    A .= 0
    @assert length(d) <= maximum(size(A))
    @inbounds @simd for i = 1:length(d)
        A[i, i] = d[i]
    end
end

function set_to_diagonal!(A::AbstractArray, d::T) where {T}
    A .= 0
    @inbounds @simd for i = 1:size(A, 1)
        A[i, i] = d
    end
end

function drop_converged_indices!(data::LOBPCGIterable{T}, residual_norms) where T
    active_indices = (residual_norms .> data.tol)
    active_indices[1:length(data.active_indices)] .&= data.active_indices
    # active_indices[sortperm(residual_norms)[1:10]].=true #[max(end - 9, 1):end]] .= true
    data.active_indices = active_indices
    
    if data.iteration <= 1
        data.active_indices .= true
    end

    # Returns the following without memory allocation
    # R = data.R[:, data.active_indices]; AR = similar(R) 
    # P = data.P[:, data.active_indices .& !data.is_restarted]
    # AP = data.AP[:, data.active_indices .& !data.is_restarted]

    R = copy_and_view(data.R, data._memoryR, data.active_indices)
    AR = view(data._memoryAR, :, 1:sum(data.active_indices))
    P = copy_and_view(data.P, data._memoryP, data.active_indices .& !data.is_restarted)
    AP = copy_and_view(data.AP, data._memoryAP, data.active_indices .& !data.is_restarted)

    return R, P, AR, AP, residual_norms
end

function copy_and_view(source, dest, active_indices)
    @assert size(source, 1) == size(dest, 1)
    @assert length(active_indices) <= min(size(source, 2), size(dest, 2)) || sum(active_indices) == 0
    indices = findall(active_indices)
    @inbounds for (counter, idx) in enumerate(indices)
        copyto!(view(dest, :, counter), view(source, :, idx))
    end
    return view(dest, :, 1:length(indices))
end

function RayleighRitz!(X::AbstractMatrix{T}, AX::AbstractMatrix{T}) where T
    XAX = Symmetric(X' * AX)
    XX = BLAS.syrk('U', 'T', one(T), X)
    λ, V = eigen!(Symmetric(X' * AX), Symmetric(XX))
    copyto!(X, X * V)
    copyto!(AX, AX * V)
    return λ, V
end

function orthonormalize!(A::AbstractMatrix{T}) where {T}
    m = size(A, 2)
    # TODO: elliminate memory allocation in the following line
    C = BLAS.syrk('U', 'T', one(T), A)
    F = cholesky!(Symmetric(C), Val(false); check=false)
    if issuccess(F)
        rdiv!(A, F.U)
        return F, true
    else
        return F, false
    end
end

function column_norms(A::AbstractMatrix{T}) where T
    v = zeros(T, size(A, 2))
	@inbounds @simd for i = 1:size(A, 2)
		v[i] = norm(view(A, :, i))
	end
	return v
end