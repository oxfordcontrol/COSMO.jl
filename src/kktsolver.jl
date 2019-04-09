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

function _kktutils_check_dims(P, A, sigma, rho)

    n = size(P, 1)
    m = size(A, 1)

    size(A,2) == n || throw(DimensionMismatch())

    length(rho)   == m || length(rho)   == 1 || throw(DimensionMismatch())
    length(sigma) == n || length(sigma) == 1 || throw(DimensionMismatch())

    return m, n
end

function _kktutils_make_kkt(P, A, sigma, rho, shape::Symbol=:F)

    S = sigma * I
    n = size(P, 1)
    m = size(A, 1)

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


# -------------------------------------
# Pardiso
# -------------------------------------

abstract type AbstractPardisoKKTSolver <: AbstractKKTSolver end

#---------------------------
#Pardiso 5.0 Direct Solver Configuration
#---------------------------

struct PardisoDirectKKTSolver <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    rhoidx::AbstractArray
    work::Vector  #working memory for Pardiso

    function PardisoDirectKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC, sigma, rho)

        m, n, K, ps, work = _pardiso_common_init(P, A, sigma, rho, PardisoSolver)
        rhoidx = _kktutils_rho_index(K, m, n)

        set_solver!(ps,Pardiso.DIRECT_SOLVER)
        _kktutils_direct_solver_init(ps, K, work, m, n)

        return new(ps, K, m, n, rhoidx, work)
    end
end

#---------------------------
#Pardiso 5.0 Indirect Solver Configuration
#---------------------------

struct PardisoIndirectKKTSolver <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    rhoidx::AbstractArray
    work::Vector  #working memory for Pardiso

    function PardisoIndirectKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC, sigma, rho)

        m, n, K, ps, work = _pardiso_common_init(P, A, sigma, rho, PardisoSolver)
        rhoidx = _kktutils_rho_index(K, m, n)

        set_solver!(ps, Pardiso.ITERATIVE_SOLVER)
        pardisoinit(ps)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        return new(ps, K, m, n, rhoidx, work)
    end
end

#---------------------------
#MKL Pardiso Solver Configuration
#---------------------------

struct MKLPardisoKKTSolver <: AbstractPardisoKKTSolver

    ps::MKLPardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    rhoidx::AbstractArray
    work::Vector  #working memory for Pardiso

    function MKLPardisoKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC, sigma, rho)

        m, n, K, ps, work = _pardiso_common_init(P, A, sigma, rho, MKLPardisoSolver)
        rhoidx = _kktutils_rho_index(K, m, n)

        _kktutils_direct_solver_init(ps, K, work, m, n)

        return new(ps ,K, m, n, rhoidx, work)
    end
end

#-------------------
#Pardiso common utils

function _pardiso_common_init(P, A, sigma, rho, Solver::Type)

    #allow Pardiso solvers do this first

    m,n  = _kktutils_check_dims(P, A, sigma, rho)
    K    = _kktutils_make_kkt(P, A, sigma, rho, :L)
    work = zeros(eltype(A), m + n)
    ps   = Solver()

    #set to symmetric indefinite
    set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)

    return m, n, K, ps, work
end


function _kktutils_direct_solver_init(ps, K, work, m, n)


    # set to allow non-default iparams
    set_iparm!(ps, 1, 1)
    # set the transpose flag (Pardiso: CSR, Julia CSC)
    set_iparm!(ps, 12, 1)
    # print output for debugging purposes
    # set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)

    # Analyze the matrix and compute a symbolic factorization.
    set_phase!(ps, Pardiso.ANALYSIS)
    pardiso(ps, K, work)

    # Compute the numeric factorization.
     set_phase!(ps, Pardiso.NUM_FACT)
     pardiso(ps, K, work)
    _pardiso_check_inertia(ps, m, n)

    ## set phase to solving for iterations
    set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
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
        _pardiso_check_inertia(s.ps, s.m, s.n)

        ## set back to solving phase
        set_phase!(s.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    end
end

function _pardiso_check_inertia(ps, m, n)

    num_pos_eigenvalues = get_iparm(ps, 22)
    num_neg_eigenvalues = get_iparm(ps, 23)

    num_neg_eigenvalues == m || throw("P matrix appears to have negative eigenvalues")
end


#number of solver threads active
get_number_of_threads(s::AbstractPardisoKKTSolver) = get_nprocs(s.ps)

#MKL can set number of threads directly
set_number_of_threads(s::MKLPardisoKKTSolver, i) = set_nprocs!(s.ps, i)

#Non-MKL can only set number of threads via an environment variable
set_number_of_threads(s::AbstractPardisoKKTSolver, i) =
    error("Pardiso 5.0 (non-MKL) must set thread count via environment variable, e.g. ENV[\"OMP_NUM_THREADS\"] = 4")
