mutable struct PsdConeTriangleLOBPCG{T} <: AbstractConvexCone{T}
    dim::Int
    sqrt_dim::Int
    iteration::Int
    is_subspace_positive::Bool
    X::Matrix{T}
    U::Matrix{T}
    lobpcg::LOBPCGIterable{T}
    exact_projections::Int
    lobpcg_iterations::Int
    max_iter::Int
    subspace_dim_history::Vector{Int}

    function PsdConeTriangleLOBPCG{T}(dim::Int; buffer_size::Int=-1, max_iter::Int=10, max_subspace_dim::Int=-1, verbosity::Int=0) where {T}
        dim >= 0       || throw(DomainError(dim, "dimension must be nonnegative"))
        sqrt_dim = Int(1 / 2 * (sqrt(8 * dim + 1) - 1)) # Solution of (sqrt_dim^2 + sqrt_dim)/2 = length(x) obtained by WolframAlpha
        sqrt_dim * (sqrt_dim + 1) / 2 == dim || throw(DomainError(dim, "dimension must be a square"))

        if buffer_size < 0
            buffer_size = Int(floor(min(max(sqrt_dim / 50, 3), 20)))
        end
        if max_subspace_dim < 0
            max_subspace_dim = Int(floor(sqrt_dim/4))
        end

        lobpcg = LOBPCGIterable(zeros(T, sqrt_dim, sqrt_dim), verbosity=verbosity,
            buffer_size=buffer_size, λ_barrier=T(0), max_dim=max_subspace_dim, min_full_iters=1)

        new(dim, sqrt_dim, 0, true,
            zeros(T, sqrt_dim, sqrt_dim), # X
            zeros(T, sqrt_dim, 0), # U
            lobpcg,
            0, 0, max_iter, # exact_projections, lobpcg_iterations
            zeros(Int, 0) # subspace_dim_history
        )
    end
end
PsdConeTriangleLOBPCG(args...; kwargs...) = PsdConeTriangleLOBPCG{DefaultFloat}(args...; kwargs...)

function get_tolerance(cone::PsdConeTriangleLOBPCG{T}) where T
    return T(max(sqrt(cone.sqrt_dim) / cone.iteration^(1.01) * 10, 1e-7))
end

function project_exact!(x::AbstractArray{T}, cone::PsdConeTriangleLOBPCG{T}) where {T}
    cone.exact_projections += 1
    populate_upper_triangle!(cone.X, x)
    λ, U  = eigen!(Symmetric(cone.X))
    zero_tol = 1e-9
    cone.is_subspace_positive = sum(λ .> zero_tol) < cone.sqrt_dim/2
    if cone.is_subspace_positive
        range = max(sum(λ .<= zero_tol) - cone.lobpcg.buffer_size + 1, 1):cone.sqrt_dim
    else
        range = 1:min(sum(λ .< -zero_tol) + cone.lobpcg.buffer_size, cone.sqrt_dim)
        λ .= -λ
        populate_upper_triangle!(cone.X, x) # Restore cone.X altered by eigen!
    end
    cone.U = view(U, :, range)
    construct_projection!(x, cone, view(λ, range))
end

function construct_projection!(x, cone::PsdConeTriangleLOBPCG{T}, λ) where {T}
    append!(cone.subspace_dim_history, size(cone.U, 2))
    scaled_U = rmul!(copy(cone.U), Diagonal(sqrt.(max.(λ, 0))))
    if cone.is_subspace_positive
        BLAS.syrk!('U', 'N', one(T), scaled_U, zero(T), cone.X)
    else
        BLAS.syrk!('U', 'N', one(T), scaled_U, one(T), cone.X)
    end
    extract_upper_triangle!(cone.X, x)
end

function project!(x::AbstractArray, cone::PsdConeTriangleLOBPCG{T}) where {T}
    sqrt_dim = cone.sqrt_dim
    cone.iteration += 1

    if size(cone.U, 2) <= cone.lobpcg.max_dim && cone.iteration > 1
        populate_upper_triangle!(cone.X, x)
        initialize!(cone.lobpcg, Symmetric(cone.X), cone.U, is_orthonormal=true,
            which = cone.is_subspace_positive ? :largest : :smallest)
        cone.lobpcg.tol = get_tolerance(cone)
        #try
            cone.lobpcg, status = lobpcg!(cone.lobpcg, cone.max_iter)
            cone.lobpcg_iterations += cone.lobpcg.iteration - 1
            #=
        catch e
            if isa(e, PosDefException)
                status = :error
                @warn "LOBPCG failed (PosDefException); switching to exact eigendecompositon."
            else
                throw(e)
            end
        end
        =#
    else
        status = :not_called
    end

    if status != :converged
        return project_exact!(x, cone)
    else
        cone.U = cone.lobpcg.X
        if !cone.is_subspace_positive
            cone.lobpcg.λ .= -cone.lobpcg.λ
        end
        construct_projection!(x, cone, cone.lobpcg.λ)
    end
end

function reset_iteration_counters!(cone::PsdConeTriangleLOBPCG)
    cone.exact_projections = 0
    cone.lobpcg_iterations = 0
    cone.iteration = 1
end

allocate_memory!(cone::PsdConeTriangleLOBPCG) = nothing