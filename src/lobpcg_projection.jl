using LinearMaps, LinearAlgebra, Random

# ----------------------------------------------------
# Positive Semidefinite Cone
# ----------------------------------------------------
mutable struct PsdConeTriangleLanczos{T} <: AbstractConvexCone{T}
    dim::Int
    n::Int
    iter_number::Int
    positive_subspace::Bool
    X::Matrix{T}
    U::Matrix{T}
    data::LOBPCGIterable{T}
    # History
    residual_history::Vector{T}
    λ_rem_history::Vector{T}
    subspace_dim_history::Vector{Int}
    λ_rem_multiplications::Vector{Int}

    function PsdConeTriangleLanczos{T}(dim::Int) where {T}
        dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
        n = Int(1 / 2 * (sqrt(8 * dim + 1) - 1)) # Solution of (n^2 + n)/2 = length(x) obtained by WolframAlpha
        n * (n + 1) / 2 == dim || throw(DomainError(dim, "dimension must be a square"))
        initial_dim = max(10, Int(floor(n / 50)))
        rng = MersenneTwister(1)
        data = LOBPCGIterable(zeros(T, n, n), verbosity = 1)
        new(dim, n, 0, true,
            zeros(T, n, n), # X
            randn(rng, T, n, initial_dim), # U
            data,
            zeros(T, 0), # residual_history
            zeros(T, 0), # λ_rem_history
            zeros(Int, 0), # subspace_dim_history
            zeros(Int, 0)
        )
    end
end
PsdConeTriangleLanczos(dim) = PsdConeTriangleLanczos{DefaultFloat}(dim)

function get_tolerance(cone::PsdConeTriangleLanczos{T}) where T
    return T(sqrt(cone.n) / cone.iter_number^(1.05)) * 10
end

function project!(x::AbstractArray, cone::PsdConeTriangleLanczos{T}) where {T}
    n = cone.n
    cone.iter_number += 1
    
    if size(cone.U, 2) < cone.data.m_max && cone.n > 50
        populate_upper_triangle!(cone.X, x)
        tol = get_tolerance(cone)
        initialize!(cone.data, cone.X, cone.U)
        cone.data.largest = cone.positive_subspace
        cone.data.tol = tol
        cone.data, status = lobpcg!(cone.data, 10)
    else
        status = :excessive_columns
    end

    if status != :converged
        return project_exact!(x, cone)
    else
        if cone.positive_subspace
            first_idx = max(sum(cone.data.λ .< 1e-6) - cone.data.buffer_size, 1)
            cone.U = cone.data.X[:, first_idx:end]
            λ = cone.data.λ[first_idx:end]
        else
            last_idx = min(sum(cone.data.λ .< -1e-6) + cone.data.buffer_size, length(cone.data.λ))
            cone.U = cone.data.X[:, 1:last_idx]
            λ = -cone.data.λ[1:last_idx]
        end
        append!(cone.subspace_dim_history, size(cone.U, 2))
        # Reconstruct projection
        scaled_U = rmul!(copy(cone.U), Diagonal(sqrt.(max.(λ, 0))))
        if cone.positive_subspace
            BLAS.syrk!('U', 'N', one(T), scaled_U, zero(T), cone.X)
        else
            BLAS.syrk!('U', 'N', one(T), scaled_U, one(T), cone.X)
        end
        extract_upper_triangle!(cone.X, x)
    end
end

function project_exact!(x::AbstractArray{T}, cone::PsdConeTriangleLanczos{T}) where {T}
    # @show "exact"
    n = cone.n

    # handle 1D case
    if length(x) == 1
        x = max.(x, zero(T))
    else
        # symmetrized square view of x
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
        if sum(λ .> 1e-6) <= sum(λ .< -1e-6)
            cone.positive_subspace = true
        else
            λ = -λ # Equivalent to considering -cone.X instead of cone.X
            cone.positive_subspace = false
        end
        idx = sum(λ .<= 1e-6) # First positive index
        # Take also a few vectors from the discarted eigenspace
        idx = max(idx - cone.data.buffer_size, 1)
        cone.U = U[:, idx:end]
    end
    append!(cone.subspace_dim_history, size(cone.U, 2))
    return nothing
end