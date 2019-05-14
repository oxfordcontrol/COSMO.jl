export QdldlKKTSolver, CholmodKKTSolver
# -------------------------------------
# abstract type defs
# -------------------------------------
abstract type AbstractKKTSolver end

# NB: all concrete examples of this type should
# implement refactor! and solve! methods and should
# implement a constructor taking (P,A,sigma,rho) as
# arguments

# -------------------------------------
# some internal utility functions
# -------------------------------------

function _kktutils_check_dims(P, A, sigma, rho)

    n = size(P, 1)
    m = size(A, 1)

    size(A,2) == n || throw(DimensionMismatch())

    length(rho)   == m || length(rho)   == 1 || throw(DimensionMismatch())
    length(sigma) == n || length(sigma) == 1 || throw(DimensionMismatch())

    return m, n
end

function _kktutils_make_kkt(P, A, sigma, rho, shape::Symbol=:F)

    S = length(sigma) == 1 ? (sigma[1]) * I : Diagonal(sigma)
    n = size(P, 1)
    m = size(A, 1)

    if length(rho) == 1
        rho = rho .* ones(m)
    end

    if     shape == :F
        #compute the full KKT matrix
        K = [P + S SparseMatrixCSC(A'); A -I]

    elseif shape == :U
        #upper triangular
        K = [triu(P) + S  SparseMatrixCSC(A'); spzeros(eltype(A), m, n)  -I]

    elseif shape == :L
        #lower triangular
        K = [tril(P)+S  spzeros(eltype(A), n, m); A  -I]

    else
        error("Bad matrix shape description")
    end

    @inbounds @simd for i = (n + 1):(n + m)
        K[i, i] = -1.0 / rho[i - n]
    end

    return K

end

function update_kkt_matrix(K, n::Int, m::Int, rho)
    if length(rho) == m
        @inbounds @simd for i = (n + 1):(n + m)
            K[i, i] = -1.0 / rho[i - n]
        end
    elseif length(rho) == 1
        @inbounds @simd for i = (n + 1):(n + m)
            K[i, i] = -1.0 / rho[1]
        end
    else
        throw(DimensionMismatch)
    end
end

# -------------------------------------
# QDLDL solver
# -------------------------------------

struct QdldlKKTSolver{Tv, Ti} <: AbstractKKTSolver

    m::Integer
    n::Integer
    ldlfact::QDLDL.QDLDLFactorisation{Tv, Ti}

    function QdldlKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma, rho) where {Tv, Ti}

        m,n = _kktutils_check_dims(P, A, sigma, rho)

        #NB: qdldl uses triu internally, but it reorders
        #with AMD first.  This way is memory inefficient
        #but saves having to permute the rhs/lhs each
        #time we solve.
        K   = _kktutils_make_kkt(P, A, sigma, rho, :F)
        ldlfact = qdldl(K)

        #check for exactly n positive eigenvalues
        positive_inertia(ldlfact) == n || error("Objective function is not convex.")

        new{Tv, Ti}(m, n, ldlfact)
    end
end

function solve!(s::QdldlKKTSolver, lhs, rhs)
    lhs .= rhs
    QDLDL.solve!(s.ldlfact, lhs)
end


function update_rho!(s::QdldlKKTSolver, rho)
    QDLDL.update_diagonal!(s.ldlfact, (s.n+1):(s.n+s.m),(-1 ./ rho))
end

# -------------------------------------
# Julia Native solver (CHOLMOD based)
# -------------------------------------
mutable struct CholmodKKTSolver{Tv, Ti} <: AbstractKKTSolver

    fact::SuiteSparse.CHOLMOD.Factor{Tv}
    K::SparseMatrixCSC{Tv, Ti}
    m::Integer
    n::Integer

    function CholmodKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma, rho) where {Tv, Ti}

         m,n  = _kktutils_check_dims(P, A, sigma, rho)
         K    = _kktutils_make_kkt(P, A, sigma, rho, :F)
         fact = ldlt(K)

        return new{Tv, Ti}(fact, K, m, n)

    end
end

solve!(s::CholmodKKTSolver, lhs, rhs) = lhs .= s.fact \ rhs

function update_rho!(s::CholmodKKTSolver, rho)
    update_kkt_matrix(s.K, s.n, s.m, rho)
    #complete restart
    s.fact = ldlt(s.K)
end

free_memory!(s::AbstractKKTSolver) = nothing
