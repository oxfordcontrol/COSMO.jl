using QDLDL, SparseArrays, Pardiso

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

function _kktutils_check_dims(P,A,sigma,rho)

    n = size(P,1)
    m = size(A,1)

    size(A,2) == n || throw(DimensionMismatch())

    length(rho)   == m || length(rho)   == 1 || throw(DimensionMismatch())
    length(sigma) == n || length(sigma) == 1 || throw(DimensionMismatch())

    return m, n
end

function _kktutils_make_kkt(P,A,sigma,rho,shape::Symbol=:F)

    R = length(rho)   == 1 ? ((1.)./rho[1])*I : Diagonal((1.)./rho)
    S = length(sigma) == 1 ? (sigma[1])*I : Diagonal(sigma)

    n = size(P,1)
    m = size(A,1)

    if     shape == :F
        #compute the full KKT matrix
        K = [P+S A'; A -R]

    elseif shape == :U
        #upper triangular
        K = [triu(P)+S  A'; spzeros(eltype(A),m,n)  -R]

    elseif shape == :L
        #lower triangular
        K = [tril(P)+S  spzeros(eltype(A),n,m); A  -R]

    else
        error("Bad matrix shape description")
    end

    return K

end

#an index into a KKT matrix indicating where
#values of 1/rho should be placed
_kktutils_rho_index(K,m,n) = diagind(K,0)[(n+1):(m+n)]




# -------------------------------------
# QDLDL solver
# -------------------------------------

struct QdldlKKTSolver <: AbstractKKTSolver

    m::Integer
    n::Integer
    ldlfact::QDLDL.QDLDLFactorisation

    function QdldlKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m,n = _kktutils_check_dims(P,A,sigma,rho)

        #NB: qdldl uses triu internally, but it reorders
        #with AMD first.  This way is memory inefficient
        #but saves having to permute the rhs/lhs each
        #time we solve.
        K   = _kktutils_make_kkt(P,A,sigma,rho,:F)
        ldlfact = qdldl(K)

        #check for exactly n positive eigenvalues
        positive_inertia(ldlfact) == n || error("Objective function is not convex.")

        new(m,n,ldlfact)
    end
end

function solve!(s::QdldlKKTSolver, lhs, rhs)
    lhs .= rhs
    QDLDL.solve!(s.ldlfact,lhs)
end


function update_rho!(s::QdldlKKTSolver, rho)
    QDLDL.update_diagonal!(s.ldlfact,(s.n+1):(s.n+s.m),(-1 ./ rho))
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

    function CholmodKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m,n  = _kktutils_check_dims(P,A,sigma,rho)
        K    = _kktutils_make_kkt(P,A,sigma,rho,:F)
        fact = ldlt(K)
        rhoidx = _kktutils_rho_index(K,m,n)

        return new(fact,K,m,n,rhoidx)

    end
end

solve!(s::CholmodKKTSolver, lhs, rhs) = lhs .= s.fact\rhs

function update_rho!(s::CholmodKKTSolver, rho)

    s.K[s.rhoidx] .= -1. ./ rho
    #complete restart
    s.fact = ldlt(s.K)

end


# -------------------------------------
# Pardiso
# -------------------------------------

abstract type AbstractPardisoKKTSolver end

#---------------------------
#Direct Solver Configuration
#---------------------------

mutable struct PardisoDirectKKTSolver <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    rhoidx::AbstractArray
    work::Vector  #working memory for Pardiso

    function PardisoDirectKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m, n, K, ps, work = _pardiso_common_init(P,A,sigma,rho)
        rhoidx = _kktutils_rho_index(K,m,n)

        set_solver!(ps,Pardiso.DIRECT_SOLVER)
        pardisoinit(ps)

        # Analyze the matrix and compute a symbolic factorization.
        set_phase!(ps, Pardiso.ANALYSIS)
        pardiso(ps, K, work)

        # Compute the numeric factorization.
        set_phase!(ps, Pardiso.NUM_FACT)
        pardiso(ps, K, work)
        _pardiso_check_inertia(ps,m,n)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        return new(ps,K,m,n,rhoidx,work)
    end
end

#---------------------------
#Indirect Solver Configuration
#---------------------------

mutable struct PardisoIndirectKKTSolver <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    rhoidx::AbstractArray
    work::Vector  #working memory for Pardiso

    function PardisoIndirectKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m, n, K, ps, work = _pardiso_common_init(P,A,sigma,rho)
        rhoidx = _kktutils_rho_index(K,m,n)

        set_solver!(ps,Pardiso.ITERATIVE_SOLVER)
        pardisoinit(ps)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        return new(ps,K,m,n,rhoidx,work)
    end
end

function _pardiso_common_init(P,A,sigma,rho)

    m,n  = _kktutils_check_dims(P,A,sigma,rho)
    K    = _kktutils_make_kkt(P,A,sigma,rho,:L)
    ps   = PardisoSolver()
    work = zeros(eltype(A),m+n)

    #set to symmetric indefinite
    set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)

    return m, n, K, ps, work
end

solve!(s::AbstractPardisoKKTSolver, lhs, rhs) = pardiso(s.ps, lhs, s.K, rhs)

function update_rho!(s::AbstractPardisoKKTSolver, rho::Union{Vector,AbstractFloat})

    s.K[s.rhoidx] .= -1. ./ rho

    #only refactor for the direct solver
    if s.ps.solver == Pardiso.DIRECT_SOLVER
        # Compute the numeric factorization again,
        #but skipping the analysis phase
        set_phase!(s.ps, Pardiso.NUM_FACT)
        pardiso(s.ps, s.K, s.work)
        _pardiso_check_inertia(s.ps,s.m,s.n)

        ## set back to solving phase
        set_phase!(s.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    end
end

function _pardiso_check_inertia(ps,m,n)

    num_pos_eigenvalues = get_iparm(ps, 22)
    num_neg_eigenvalues = get_iparm(ps, 23)

    num_neg_eigenvalues == m || throw("P matrix appears to have negative eigenvalues")
end
