using LinearAlgebra, Statistics
using Printf, Random

mutable struct LOBPCGIterable{T}
    n::Int
    m::Int
    m_initial::Int
    iteration::Int
    A::AbstractMatrix{T}
    λ_barrier::T
    X::Matrix{T}
    P::Matrix{T}

    AX::Matrix{T}
    AP::Matrix{T}
    residual_norms::Vector{T}

    λ::Vector{T}
    condition_estimates::Vector{T}
    active_indices::BitVector
    tol::T
    verbosity::Int
    buffer_size::Int
    delta::Int
    which::Symbol
    min_full_iters::Int
    max_dim::Int
    print_interval::Int

    start_time::Float64

    function LOBPCGIterable(A::AbstractMatrix{T};
        tol::T=T(1e-4), verbosity=1,
        buffer_size::Int=max(Int(floor(size(A, 1) / 50)), 3),
        which=:largest,
        λ_barrier=T(NaN),
        max_dim=floor(size(A, 1)/4),
        min_full_iters=0,
        print_interval=1
    ) where T
        n = size(A, 1)

        @assert issymmetric(A) "Matrix must be symmetric"
        @assert which == :largest || which == :smallest "The argument 'which' can take the values :largest or :smallest" 
        new{T}(n, 0, 0, 0, A,
            which == :largest ? λ_barrier : -λ_barrier,
            Matrix{T}(undef, n, 0), Matrix{T}(undef, n, 0), # X, P
            Matrix{T}(undef, n, 0), Matrix{T}(undef, n, 0), # ΑΧ, AP
            zeros(T, 0), # residual_norms
            zeros(T, 0), # λ
            [-log10(eps(T))/2], # condition_estimates
            BitArray(undef, 0), # active_indices
            tol, verbosity,
            buffer_size, # buffer
            0,
            which,
            min_full_iters, max_dim,
            print_interval,
            0.0
        )
    end
end

function initialize!(data::LOBPCGIterable{T}, A::AbstractMatrix{T}, X::Matrix{T}; kwargs...) where T
    @assert data.n == size(A, 1) == size(A, 2)
    @assert issymmetric(A) "Matrix must be symmetric"
    data.A = A
    initialize!(data, X; kwargs...)
end

function initialize!(data::LOBPCGIterable{T}, X::Matrix{T}; is_orthonormal=false, which=:unchanged) where T
    data.m = size(X, 2)
    data.m_initial = size(X, 2)
    @assert size(X, 1) == data.n
    if which != :unchanged
        @assert which == :largest || which == :smallest "The argument 'which' can take the values :largest or :smallest" 
        data.which = which
    end
    data.m > data.max_dim && throw(ArgumentError("Too many columns in initial matrix.")) # m_max?

    if data.m < data.buffer_size
        data.buffer_size = data.m
    end
    data.X = X
    data.AX = mul_A(data.X, data)
    restart_subspace!(data)
    data.λ, _ = rayleigh_ritz!(data.X, data.AX; is_orthonormal=is_orthonormal)
    data.active_indices = BitArray(undef, data.m)
    data.active_indices .= true

    append!(data.residual_norms, T(NaN))
    data.iteration = 0
    data.condition_estimates = [-log10(eps(T))/2]
end

function print_statistics(data::LOBPCGIterable)
    if data.iteration == 1 && data.verbosity > 0
        println("-----------------------")
        println("Computing ", string(data.which), " eigenpairs of a matrix of size: ", size(data.A))
        println("Number of initial Ritz pairs: ", size(data.X, 2), ". Tolerance: ", data.tol)
        println("Maximum allowed subspace size: ", data.max_dim, ". Minimum full iterations: ", data.min_full_iters)
        if !isnan(data.λ_barrier)
            range_string = data.which == :largest ? string("(", data.λ_barrier,", Inf)") : string("(-Inf, ", -data.λ_barrier,")")
            println("All the eigenvalues in ", range_string, " were requested. Buffer size: ", data.buffer_size)
        end
        println("-----------------------")
    end
    if data.verbosity > 1
        if mod(data.iteration - 1, data.print_interval*10) == 0
            @printf("Iter\t Non-conv pairs\t Residuals(min|max)\t Ritz Values(min|λmax) \t Time(s) \t Restarted\n")
        end
        if data.delta == 0
            size_string = string(sum(data.active_indices), "      ")
        else
            size_string = string(sum(data.active_indices), " (", data.delta, ")  ")
        end
        min_λ = data.which == :largest ? minimum(data.λ) : -maximum(data.λ)
        max_λ = data.which == :largest ? maximum(data.λ) : -minimum(data.λ)
        if mod(data.iteration - 1, data.print_interval) == 0
            @printf("%d\t %s\t %.2e | %.2e\t %+.2e | %+.2e \t %.3e \t %s\n",
                data.iteration, size_string, minimum(data.residual_norms), maximum(data.residual_norms),
                min_λ, max_λ, time() - data.start_time,
                is_restarted(data) ? "YES" : "-")
        end
    end
end

function mul_A!(Y, X, data::LOBPCGIterable)
    mul!(Y, data.A, X)
    if data.which == :smallest
        @. Y = -Y
    end
    return Y
end
mul_A(X, data::LOBPCGIterable) = return mul_A!(similar(X), X, data)

function lobpcg!(data::LOBPCGIterable{Float64}, max_iter=10)
    data.start_time = time()
    status = :not_converged
    while data.iteration < max_iter && status == :not_converged
        data.iteration += 1
        status = iterate!(data)
    end

    if data.verbosity > 0
        if status == :excessive_column_num
            println("Too many desired eigenvalues. Use direct eigendecomposition instead.")
        elseif status == :excessive_column_num
        end
        @printf("%d eigenvalues converged in %d iterations, %.3e seconds\n", sum(.!data.active_indices), data.iteration - 1, time() - data.start_time)
    end
    if data.which == :smallest
        @. data.λ = -data.λ
    end
    return data, status
end

function lobpcg(A::AbstractMatrix{Float64}, X::Matrix{Float64}, max_iter=10;  kwargs...)
    data = LOBPCGIterable(A; kwargs...)
    initialize!(data, A, X)
    lobpcg!(data, max_iter)
end

function iterate!(data::LOBPCGIterable{T}) where T
    R = data.AX - data.X*Diagonal(data.λ)
    R, data.P, data.AP, converged = drop_converged_indices(data, R)
    converged && return :converged
    R -= data.X*(data.X'*R) # Orthogonalize w.r.t. X
    _, success = orthonormalize!(R)
    if !success 
        @warn "LOBPCG: Failed to orthonormalize residual"
        return :ill_conditioned
    end

    ### Main computational time in the following time
    AR = mul_A(R, data)
    ###
    F, success = orthonormalize!(data.P)
    if success
        rdiv!(data.AP, F.U)
    else
        data.verbosity > 0 && @warn string("Orthonormalization of P on iteration ", data.iteration, " failed. Restarting subspace.")
        restart_subspace!(data)
    end
    
    # The explicit_computation flag is not consistent with Knyzev's MATLAB's implementation.
    # As a matter of fact, Knyazev has two different MATLAB versions one of which has explicit_computation always true. See
    # https://github.com/lobpcg/blopex/blob/3e5ad076eb41aad0d7d4204bc7c70c151dedaaf9/blopex_tools/matlab/lobpcg/lobpcg_Rev_4_13.m
    # https://github.com/lobpcg/blopex/blob/33b5e9f7b9de3a9c28edbc0e22381f3f2edce315/blopex_matlab/driver/lobpcg.m
    explicit_computation = any(data.residual_norms .<= eps(T)^T(0.5))
    G, Q = compute_gramm_matrix(data, R, AR, explicit_computation)
    success, condition_estimate = check_conditioning(data, Q, R)
    if !success && !is_restarted(data)
        data.verbosity > 1 && @info string("Restarting due to condition estimate: ", condition_estimate)
        G, Q = restart_subspace!(data, G, Q)
        success, condition_estimate = check_conditioning(data, Q, R)
    end
    !success && data.verbosity > 0 && @warn "Gramm matrix ill-conditioned: results unpredictable"
    push!(data.condition_estimates, condition_estimate)

    ### Main computational time in the following block
    if !explicit_computation && is_restarted(data)
        λ, W = eigen!(G) # Q is simply identity in this case
    else
        λ, W = eigen!(G, Q)
    end
    ###
    print_statistics(data)
    return update_subspaces!(data, R, AR, λ, W)
end

function check_conditioning(data, Q, R)
    condition_estimate = log10(cond(Q)) + 1;
    condition_estimate_mean = mean(data.condition_estimates[
        max(1, data.iteration - 10 - Int(round(log(size(R, 2))))):data.iteration])
    failure = (condition_estimate / condition_estimate_mean > 2 && condition_estimate > 2) || condition_estimate > 8
    return !failure, condition_estimate
end

is_restarted(data) = length(data.P) == 0

function restart_subspace!(data::LOBPCGIterable{T}) where T
    data.P = zeros(T, data.n, 0)
    data.AP = zeros(T, data.n, 0)
end

function restart_subspace!(data::LOBPCGIterable{T}, G, Q) where T
    new_dim = size(G, 1) - size(data.P, 2)
    restart_subspace!(data)
    return Symmetric(view(G.data, 1:new_dim, 1:new_dim)), Symmetric(view(Q.data, 1:new_dim, 1:new_dim))
end

function update_subspaces!(data::LOBPCGIterable{T}, R, AR, λ, W) where T
    subspace_diff = count(λi -> λi > data.λ_barrier, λ) - data.m

    # Contract subspace if necessary
    data.delta = 0
    if subspace_diff < -data.buffer_size && !isnan(data.λ_barrier)
        m_contracted = max(data.m + subspace_diff, data.buffer_size)
        data.delta = m_contracted - data.m
        data.m = m_contracted
        data.active_indices = data.active_indices[-data.delta+1:end]
    end
    data.λ = λ[end-data.m+1:end]
    W1 = view(W, :, size(W, 2) - data.m + 1:size(W,2))

    W11 = view(W1, 1:size(data.X, 2), :); W21 = view(W1, size(data.X, 2)+1:size(W, 1), :)
    data.P = [R data.P] * W21
    data.X = data.X * W11 + data.P
    data.AP = [AR data.AP] * W21
    data.AX = data.AX * W11 + data.AP

    # Expand subspace if necessary
    if subspace_diff > -data.buffer_size/2 && !isnan(data.λ_barrier) 
        data.delta = subspace_diff + data.buffer_size
        data.m += data.delta
        data.m > data.max_dim && return :excessive_column_num
        X_ = randn(data.n, data.m - size(data.X, 2))
        data.X = [X_ data.X]
        # @show svdvals(data.X)
        data.AX = [mul_A(X_, data) data.AX]
        data.λ, _ = rayleigh_ritz!(data.X, data.AX)
        # @show norm(data.X'*data.X - I)
        prepend!(data.active_indices, [true for i = 1:data.delta])
    end

    return :not_converged
end

function compute_gramm_matrix(data::LOBPCGIterable{T}, R, AR, explicit_computation) where T
    #= Computes the following two matrices efficiently
    G = Symmetric([data.X R data.P]' * mul_A([data.X R data.P], data))
    Q = Symmetric([data.X R data.P]' * [data.X R data.P])
    =#
    X, P, AX, AP = data.X, data.P, data.AX, data.AP
    if !explicit_computation
        XAX = diagm(0 => data.λ)
        XX, RR, PP = I, I, I
        XR = zeros(T, size(X, 2), size(R, 2))
    else
        XAX = X'*AX
        XX, RR, PP = X'*X, R'*R, P'*P
        XR = X'*R
    end
    XAR = X'*AR;    XAP = X'*AP
    RAR = R'*AR;    RAP = R'*AP
                    PAP = P'*AP
    XP = X'*P; RP = R'*P

    G = [XAX    XAR     XAP;
         XAR'   RAR     RAP;
         XAP'   RAP'    PAP]

    Q = [XX     XR      XP;
         XR'    RR      RP;
         XP'    RP'     PP]
    #@assert norm(Symmetric(Q) - Symmetric([data.X R data.P]' * [data.X R data.P])) <= 1e-6 norm(Symmetric(Q) - Symmetric([data.X R data.P]' * [data.X R data.P]))
    return Symmetric(G), Symmetric(Q)
end

function drop_converged_indices(data::LOBPCGIterable{T}, R) where T
    data.residual_norms = column_norms(R)
    # Technically, the &= in the next command should be =
    # We have it like this to match the original implementation in MATLAB (BLOPEX).
    @. data.active_indices &= (data.residual_norms > data.tol)
    converged = false
    if data.iteration > data.min_full_iters
        R = R[:, data.active_indices]
        if !is_restarted(data)
            data.P = data.P[:, data.active_indices[1:size(data.P, 2)]]
            data.AP = data.AP[:, data.active_indices[1:size(data.AP, 2)]]
        end
        num_converged = length(data.active_indices) - count(data.active_indices)
        if num_converged >= data.m_initial
            # Terminate if all eigeivnvalues in the desired region have converged
            converged = !any(data.active_indices[data.λ .+ data.residual_norms .> data.λ_barrier])
        end
        converged |= num_converged == data.m
    end
    
    return R, data.P, data.AP, converged
end

function rayleigh_ritz!(X::AbstractMatrix{T}, AX::AbstractMatrix{T}; is_orthonormal=false) where T
    if is_orthonormal
        λ, V = eigen!(Symmetric(X' * AX))
    else
        λ, V = eigen!(Symmetric(X' * AX), Symmetric(X'*X))
    end
    copyto!(X, X * V)
    copyto!(AX, AX * V)
    return λ, V
end

function orthonormalize!(A::AbstractMatrix{T}) where {T}
    # @show svdvals(A)
    m = size(A, 2)
    F = cholesky!(Symmetric(A'*A), Val(false); check=false)
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