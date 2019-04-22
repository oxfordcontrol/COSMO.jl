using LinearMaps, IterativeSolvers

mutable struct CGKKTSolver{TP, TA, T} <: COSMO.AbstractKKTSolver
    m::Integer
    n::Integer
    P::TP
    A::TA
    σ::T
    ρ::Vector{T}
    # Memory and statistics for cg
    previous_cg_solution::Vector{T}
    counter::Int
    multiplications::Vector{Int}
    tmp_n::Vector{T}
    tmp_m::Vector{T}

    function CGKKTSolver(P::TP, A::TA, σ::T, ρ) where {TP, TA, T}
        m, n = size(A)
        new{TP, TA, T}(m, n, P, A, σ, ρ,
            zeros(T, n), 1, zeros(Int, 1), zeros(T, n), zeros(T, m))
    end
end

function solve!(S::CGKKTSolver, y::AbstractVector{T}, x::AbstractVector{T}) where {T}
    # Solves the (KKT) linear system
    # | P + σI     A'  |  |y1|  =  |x1|
    # | A        -I/ρ  |  |y2|  =  |x2|
    # x1,y1: R^n, x2/y2: R^m
    # where [y1; y2] := y, [x1; x2] := x 

    # In particular we perform cg on the reduced system
    # (P + σΙ + Α'ρΑ)y1 = x1 + A'ρx2
    # And then recover y2 as
    # y2 = ρ(Ay1 - x2)

    x1 = view(x, 1:S.n); y1 = view(y, 1:S.n)
    x2 = view(x, S.n+1:S.n+S.m); y2 = view(y, S.n+1:S.n+S.m)

    # Form right-hand side for cg: y1 = x1 + A'ρx2
    @. y2 = S.ρ*x2
    mul!(y1, S.A', y2)
    y1 .+= x1
    function cg_mul!(y::AbstractVector{T}, x::AbstractVector{T}) where {T}
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
    L = LinearMap(cg_mul!, S.n; ismutating=true, issymmetric=true)
    tol = 1/S.counter^2
    cg!(S.previous_cg_solution, L, y1, tol=tol/norm(y1), initially_zero=false) #, maxiter=10)
    # @assert tol > norm(L*S.previous_cg_solution - y1)
    copyto!(y1, S.previous_cg_solution)

    # y2 = Ay1 - x2
    mul!(y2, S.A, y1)
    axpy!(-one(T), x2, y2)
    @. y2 .*= S.ρ

    S.counter += 1
    push!(S.multiplications, 0)

    return y
end

function update_rho!(S::CGKKTSolver, ρ)
    S.ρ .= ρ
end