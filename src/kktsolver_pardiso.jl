using Pardiso: Pardiso, set_solver!, pardisoinit, set_phase!, set_iparm!, set_msglvl!, get_iparm, pardiso, PardisoSolver, MKLPardisoSolver, get_nprocs, set_nprocs!, set_matrixtype!
export PardisoDirectKKTSolver, PardisoIndirectKKTSolver, MKLPardisoSolver

abstract type AbstractPardisoKKTSolver <: AbstractKKTSolver end

#---------------------------
#Pardiso 5.0 Direct Solver Configuration
#---------------------------

struct PardisoDirectKKTSolver{Tv, Ti} <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC{Tv, Ti}
    m::Integer
    n::Integer
    work::Vector{Tv}  #working memory for Pardiso

    function PardisoDirectKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma, rho;
        iparm::Dict{Int64, Int64} = Dict{Int64, Int64}(), msg_level_on::Bool = false) where {Tv, Ti}

        m, n, K, ps, work = _pardiso_common_init(P, A, sigma, rho, PardisoSolver, iparm, msg_level_on)

        set_solver!(ps,Pardiso.DIRECT_SOLVER)
        _kktutils_direct_solver_init(ps, K, work, m, n)

        return new{Tv, Ti}(ps, K, m, n, work)
    end
end

#---------------------------
#Pardiso 5.0 Indirect Solver Configuration
#---------------------------

struct PardisoIndirectKKTSolver{Tv, Ti} <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC{Tv, Ti}
    m::Integer
    n::Integer
    work::Vector{Tv}  #working memory for Pardiso

    function PardisoIndirectKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma, rho;
        iparm::Dict{Int64, Int64} = Dict{Int64, Int64}(), msg_level_on::Bool = false) where {Tv, Ti}

        m, n, K, ps, work = _pardiso_common_init(P, A, sigma, rho, PardisoSolver, iparm, msg_level_on)

        set_solver!(ps, Pardiso.ITERATIVE_SOLVER)
        pardisoinit(ps)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        return new{Tv, Ti}(ps, K, m, n, work)
    end
end

#---------------------------
#MKL Pardiso Solver Configuration
#---------------------------

struct MKLPardisoKKTSolver{Tv, Ti} <: AbstractPardisoKKTSolver

    ps::MKLPardisoSolver
    K::SparseMatrixCSC{Tv, Ti}
    m::Integer
    n::Integer
    work::Vector{Tv}  #working memory for Pardiso

    function MKLPardisoKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma, rho;
        iparm::Dict{Int64, Int64} = Dict{Int64, Int64}(), msg_level_on::Bool = false, num_threads::Int64 = 1) where {Tv, Ti}

        m, n, K, ps, work = _pardiso_common_init(P, A, sigma, rho, MKLPardisoSolver, iparm, msg_level_on)

        _kktutils_direct_solver_init(ps, K, work, m, n)
        if num_threads >= 1
            set_nprocs!(ps, num_threads)
        else
            error("The number of CPU threads has to be greater than 0.")
        end

        return new{Tv, Ti}(ps ,K, m, n, work)
    end
end

#-------------------
#Pardiso common utils

function _pardiso_common_init(P, A, sigma, rho, Solver::Type, iparm::Dict{Int64, Int64}, msg_level_on::Bool)

    #allow Pardiso solvers do this first

    m,n  = _kktutils_check_dims(P, A, sigma, rho)
    K    = _kktutils_make_kkt(P, A, sigma, rho, :L)
    work = zeros(eltype(A), m + n)
    ps   = Solver()

    #set to symmetric indefinite
    set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)

    if msg_level_on
        set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    else
        set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_OFF)
    end

    # set to allow non-default iparams
    set_iparm!(ps, 1, 1)

    # set user specific parameters
    for (i, v) in iparm
        set_iparm!(ps, i, v)
    end


    return m, n, K, ps, work
end


function _kktutils_direct_solver_init(ps, K, work, m, n)


    # set the transpose flag (Pardiso: CSR, Julia CSC)
    set_iparm!(ps, 12, 1)

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
    update_kkt_matrix(s.K, s.n, s.m, rho)

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

function free_memory!(s::AbstractPardisoKKTSolver)
    set_phase!(s.ps, Pardiso.RELEASE_ALL)
    pardiso(s.ps)
end

#number of solver threads active
get_number_of_threads(s::AbstractPardisoKKTSolver) = get_nprocs(s.ps)

#MKL can set number of threads directly
set_number_of_threads(s::MKLPardisoKKTSolver, i) = set_nprocs!(s.ps, i)

#Non-MKL can only set number of threads via an environment variable
set_number_of_threads(s::AbstractPardisoKKTSolver, i) =
    error("Pardiso 5.0 (non-MKL) must set thread count via environment variable, e.g. ENV[\"OMP_NUM_THREADS\"] = 4")
