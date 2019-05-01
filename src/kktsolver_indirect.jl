using LinearMaps, IterativeSolvers

mutable struct IndirectReducedKKTSolver{TP, TA, T} <: COSMO.AbstractKKTSolver
    m::Integer
    n::Integer
    P::TP
    A::TA
    σ::T
    ρ::Vector{T}
    # Memory and statistics for cg
    solver::Symbol
    previous_cg_solution::Vector{T}
    iteration_counter::Int
    multiplications::Vector{Int}
    tmp_n::Vector{T}
    tmp_m::Vector{T}

    function IndirectReducedKKTSolver(P::TP, A::TA, σ::T, ρ; solver::Symbol=:CG) where {TP, TA, T}
        m, n = size(A)
        if isa(ρ, T)
            ρ = ρ*ones(T, m)
        end
        @assert eltype(P) == eltype(A) == T "Inconsistent element types."
        @assert solver == :CG || solver == :MINRES "Solver symbol must be either :CG or :MINRES"
        new{TP, TA, T}(m, n, P, A, σ, ρ,
            solver,
            zeros(T, n), 1, zeros(Int, 1), zeros(T, n), zeros(T, m))
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

    # Form right-hand side for cg: y1 = x1 + A'ρx2
    @. y2 = S.ρ*x2
    mul!(y1, S.A', y2)
    y1 .+= x1
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
    if S.solver == :CG
        cg!(S.previous_cg_solution, L, y1, tol=get_tolerance(S)/norm(y1)) #, maxiter=10)
    elseif S.solver == :MINRES
        init_residual = norm(L*S.previous_cg_solution - y1)
        minres!(S.previous_cg_solution, L, y1, tol=get_tolerance(S)/init_residual) #, maxiter=10)
    end
    # @assert get_tolerance(S) > norm(L*S.previous_cg_solution - y1)
    copyto!(y1, S.previous_cg_solution)

    # y2 = Ay1 - x2
    mul!(y2, S.A, y1)
    axpy!(-one(T), x2, y2)
    @. y2 .*= S.ρ

    S.iteration_counter += 1
    push!(S.multiplications, 0)

    return y
end

function update_rho!(S::IndirectReducedKKTSolver, ρ)
    S.ρ .= ρ
end

function get_tolerance(S::IndirectReducedKKTSolver)
    return 1/S.iteration_counter^2
end