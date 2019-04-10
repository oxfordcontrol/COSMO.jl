
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

#an index into a KKT matrix indicating where
#values of 1/rho should be placed
_kktutils_rho_index(K, m, n) = diagind(K, 0)[(n+1):(m+n)]




# -------------------------------------
# QDLDL solver
# -------------------------------------

struct QdldlKKTSolver <: AbstractKKTSolver

    m::Integer
    n::Integer
    ldlfact::QDLDL.QDLDLFactorisation

    function QdldlKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC, sigma, rho)

        m,n = _kktutils_check_dims(P, A, sigma, rho)

        #NB: qdldl uses triu internally, but it reorders
        #with AMD first.  This way is memory inefficient
        #but saves having to permute the rhs/lhs each
        #time we solve.
        K   = _kktutils_make_kkt(P, A, sigma, rho, :F)
        ldlfact = qdldl(K)

        #check for exactly n positive eigenvalues
        positive_inertia(ldlfact) == n || error("Objective function is not convex.")

        new(m, n, ldlfact)
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
mutable struct CholmodKKTSolver <: AbstractKKTSolver

    fact::SuiteSparse.CHOLMOD.Factor
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    rhoidx::AbstractArray

    function CholmodKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC, sigma, rho)

         m,n  = _kktutils_check_dims(P, A, sigma, rho)
         K    = _kktutils_make_kkt(P, A, sigma, rho, :F)
         fact = ldlt(K)
         rhoidx = _kktutils_rho_index(K, m, n)

        return new(fact, K, m, n, rhoidx)

    end
end

solve!(s::CholmodKKTSolver, lhs, rhs) = lhs .= s.fact \ rhs

function update_rho!(s::CholmodKKTSolver, rho)

    s.K[s.rhoidx] .= -1. ./ rho
    #complete restart
    s.fact = ldlt(s.K)

end