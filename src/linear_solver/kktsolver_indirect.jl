using .LinearMaps, .IterativeSolvers
export CGIndirectKKTSolver, MINRESIndirectKKTSolver
mutable struct IndirectReducedKKTSolver{TP, TA, T} <: AbstractKKTSolver
    m::Integer
    n::Integer
    P::TP
    A::TA
    σ::T
    ρ::Vector{T}
    tol_constant::T
    tol_exponent::T
    # Memory and statistics for cg
    solver_type::Symbol
    previous_solution::Vector{T}
    iteration_counter::Int
    multiplications::Vector{Int}
    tmp_n::Vector{T} # n-dimensional used to avoid memory allocation
    tmp_m::Vector{T} # m-dimensional used to avoid memory allocation

    function IndirectReducedKKTSolver(P::TP, A::TA, σ::T, ρ;
        solver_type::Symbol=:CG, tol_constant::T=T(1.0), tol_exponent::T=T(1.5)
        ) where {TP, TA, T}
        m, n = size(A)
        if isa(ρ, T)
            ρ = ρ*ones(T, m)
        end
        @assert eltype(P) == eltype(A) == T "Inconsistent element types."
        @assert solver_type == :CG || solver_type == :MINRES "Solver symbol must be either :CG or :MINRES"
        new{TP, TA, T}(m, n, P, A, σ, ρ,
            tol_constant, tol_exponent,
            solver_type,
            zeros(T, n), 1, zeros(Int, 0), zeros(T, n), zeros(T, m))
    end
end

function solve!(S::IndirectReducedKKTSolver, y::AbstractVector{T}, x::AbstractVector{T}) where {T}
    # Solves the (KKT) linear system
    # | P + σI     A'  |  |y1|  =  |x1|
    # | A        -I/ρ  |  |y2|  =  |x2|
    # x1,y1: R^n, x2/y2: R^m
    # where [y1; y2] := y, [x1; x2] := x

    # In particular we perform cg/minres on the reduced system
    # (P + σΙ + Α'ρΑ)y1 = x1 + A'ρx2
    # And then recover y2 as
    # y2 = ρ(Ay1 - x2)

    x1 = view(x, 1:S.n); y1 = view(y, 1:S.n)
    x2 = view(x, S.n+1:S.n+S.m); y2 = view(y, S.n+1:S.n+S.m)

    # Form right-hand side for cg/minres: y1 = x1 + A'ρx2
    @. y2 = S.ρ*x2
    mul!(y1, S.A', y2)
    y1 .+= x1

    push!(S.multiplications, 0)
    function reduced_mul!(y::AbstractVector{T}, x::AbstractVector{T}) where {T}
        # y = (P + σΙ + Α'ρΑ)x
        mul!(S.tmp_m, S.A, x)
        @. S.tmp_m .*= S.ρ
        mul!(S.tmp_n, S.A', S.tmp_m)
        axpy!(S.σ, x, S.tmp_n)
        mul!(y, S.P, x)
        axpy!(one(T), S.tmp_n, y)
        S.multiplications[end] += 1
        return y
    end
    L = LinearMap(reduced_mul!, S.n; ismutating=true, issymmetric=true)
    if S.solver_type == :CG
        cg!(S.previous_solution, L, y1, abstol=get_tolerance(S)/norm(y1), reltol = 0.)
    elseif S.solver_type == :MINRES
        init_residual = norm(L*S.previous_solution - y1)
        minres!(S.previous_solution, L, y1, abstol=get_tolerance(S)/init_residual, reltol = 0.)
    end
    # Sanity check for tolerance
    # might not always hold for MINRES, as its termination criterion is approximate, (see its documentation)
    # @assert get_tolerance(S) > norm(L*S.previous_solution - y1)
    copyto!(y1, S.previous_solution)

    # y2 = Ay1 - x2
    mul!(y2, S.A, y1)
    axpy!(-one(T), x2, y2)
    @. y2 .*= S.ρ

    S.iteration_counter += 1

    return y
end

mutable struct IndirectKKTSolver{TP, TA, T} <: AbstractKKTSolver
    m::Integer
    n::Integer
    P::TP
    A::TA
    σ::T
    ρ::Vector{T}
    tol_constant::T
    tol_exponent::T
    # Memory and statistics for cg
    solver_type::Symbol
    previous_solution::Vector{T}
    iteration_counter::Int
    multiplications::Vector{Int}
    tmp_n::Vector{T} # n-dimensional used to avoid memory allocation
    tmp_m::Vector{T} # m-dimensional used to avoid memory allocation

    function IndirectKKTSolver(P::TP, A::TA, σ::T, ρ; solver_type::Symbol=:MINRES,
        tol_constant::T=T(1.0), tol_exponent::T=T(1.5)
        ) where {TP, TA, T}
        m, n = size(A)
        if isa(ρ, T)
            ρ = ρ*ones(T, m)
        end
        @assert eltype(P) == eltype(A) == T "Inconsistent element types."
        @assert solver_type == :MINRES "Solver symbol must be :MINRES"
        new{TP, TA, T}(m, n, P, A, σ, ρ,
            tol_constant, tol_exponent,
            solver_type,
            zeros(T, m + n), 1, zeros(Int, 0), zeros(T, n), zeros(T, m))
    end
end

function solve!(S::IndirectKKTSolver, y::AbstractVector{T}, x::AbstractVector{T}) where {T}
    # Solves the (KKT) linear system
    # | P + σI     A'  |  |y1|  =  |x1|
    # | A        -I/ρ  |  |y2|  =  |x2|
    # x1,y1: R^n, x2/y2: R^m
    # where [y1; y2] := y, [x1; x2] := x
    push!(S.multiplications, 0)
    function kkt_mul!(y::AbstractVector{T}, x::AbstractVector{T}) where {T}
        # Performs
        # |y1| = | P + σI     A'  |  |x1|
        # |y2| = | A        -I/ρ  |  |x2|
        x1 = view(x, 1:S.n); y1 = view(y, 1:S.n)
        x2 = view(x, S.n+1:S.n+S.m); y2 = view(y, S.n+1:S.n+S.m)

        mul!(S.tmp_n, S.A', x2)
        axpy!(S.σ, x1, S.tmp_n)
        mul!(y1, S.P, x1)
        axpy!(one(T), S.tmp_n, y1)

        @. y2 = -x2/S.ρ
        mul!(S.tmp_m, S.A, x1)
        axpy!(one(T), S.tmp_m, y2)

        S.multiplications[end] += 1
        return y
    end
    L = LinearMap(kkt_mul!, S.n + S.m; ismutating=true, issymmetric=true)
    if S.solver_type == :MINRES
        init_residual = norm(L*S.previous_solution - x)
        minres!(S.previous_solution, L, x, abstol=get_tolerance(S)/init_residual, reltol = 0.)
    end
    # Sanity check for tolerance
    # might not always hold for MINRES, as its termination criterion is approximate, (see its documentation)
    # @assert get_tolerance(S) > norm(L*S.previous_solution - y)
    copyto!(y, S.previous_solution)

    S.iteration_counter += 1

    return y
end

function update_rho!(S::Union{IndirectReducedKKTSolver, IndirectKKTSolver}, ρ)
    S.ρ .= ρ
end

function get_tolerance(S::Union{IndirectReducedKKTSolver, IndirectKKTSolver})
    return S.tol_constant/S.iteration_counter^S.tol_exponent
end

# For convenience we define CG and Minres wrappers around IndirectKKTSolver and IndirectReducedKKTSolver
struct CGIndirectKKTSolver{TP, TA, T} <: AbstractKKTSolver
    indirect_kktsolver::IndirectReducedKKTSolver{TP, TA, T}
    function CGIndirectKKTSolver(P::TP, A::TA, σ::T, ρ; tol_constant::T=T(1.0), tol_exponent::T=T(1.5)) where {TP, TA, T}
        new{TP, TA, T}(IndirectReducedKKTSolver(P, A, σ, ρ; solver_type = :CG, tol_constant = tol_constant, tol_exponent = tol_exponent))
    end
end

struct MINRESIndirectKKTSolver{TP, TA, T} <: AbstractKKTSolver
    indirect_kktsolver::IndirectKKTSolver{TP, TA, T}
    function MINRESIndirectKKTSolver(P::TP, A::TA, σ::T, ρ; tol_constant::T=T(1.0), tol_exponent::T=T(1.5)) where {TP, TA, T}
        new{TP, TA, T}(IndirectKKTSolver(P, A, σ, ρ; solver_type = :MINRES, tol_constant = tol_constant, tol_exponent = tol_exponent))
    end
end

update_rho!(S::Union{CGIndirectKKTSolver, MINRESIndirectKKTSolver}, ρ) = update_rho!(S.indirect_kktsolver, ρ)
get_tolerance(S::Union{CGIndirectKKTSolver, MINRESIndirectKKTSolver}, ρ) = get_tolerance(S.indirect_kkt_solver, ρ)
solve!(S::Union{CGIndirectKKTSolver, MINRESIndirectKKTSolver}, y::AbstractVector{T}, x::AbstractVector{T}) where {T} = solve!(S.indirect_kktsolver, y, x)
