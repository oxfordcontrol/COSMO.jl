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

    #short circuit for fast triu form
    if shape == :U
        K = _kkt_fast_uppertriu(P, A, sigma, rho)
        return K
    end

    n = size(P, 1)
    m = size(A, 1)
    S = length(sigma) == 1 ? (sigma[1]) * I : Diagonal(sigma)
    rhoinv  = Vector{Float64}(undef,m)
    rhoinv .= (-1.0./rho)
    D       = SparseMatrixCSC(Diagonal(rhoinv))

    if  shape == :F
        #compute the full KKT matrix
        K = [P + S SparseMatrixCSC(A'); A D]

    elseif shape == :L
        #lower triangular
        K = [tril(P)+S  spzeros(eltype(A), n, m); A  D]

    else
        error("Bad matrix shape description")
    end

    return K

end

function _kkt_fast_uppertriu(P::SparseMatrixCSC{Tf, Ti}, A::SparseMatrixCSC{Tf, Ti}, sigma, rho) where {Tf, Ti}

    n = size(P, 1)
    m = size(A, 1)

    if length(sigma) == 1
        S = Vector{typeof(sigma)}(undef,n)
        fill!(S,sigma)
    else
        S = sigma
    end

    if length(rho) == 1
        R  = Array{typeof(rho)}(undef,m)
        fill!(R,-1/rho)
    else
        R  = Array{typeof(rho[1])}(undef,m)
        R .= (-1.0./rho)
    end

    #make the column counter of K
    Kcolnz = zeros(Ti,m+n)

    #work out how many nonzeros are in K.   First
    #count nonzeros in the strict upper triangle of (P)
    @inbounds for cidx = 1:n
        for j = (P.colptr[cidx]):(P.colptr[cidx+1]-1)
            ridx = P.rowval[j]
            if (ridx < cidx)
                Kcolnz[cidx] += 1
            else
                break
            end
        end
    end
    #count the nonzeros in columns in A'
    @inbounds for idx in A.rowval
        Kcolnz[idx + n] += 1
    end

    #every element on the diagonal is also nonzero
    Kcolnz .+= 1

    #make the column pointer and nonzero count
    Kp   = Vector{Ti}(undef,m+n+1)
    Kp[1] = 1
    @inbounds for i = 1:(m+n)
        Kp[i+1] = Kp[i] + Kcolnz[i]
    end

    KNextInCol = copy(Kp);        #next value in column goes here
    nnzK = Kp[end] - 1
    Ki   = Vector{Ti}(undef,nnzK)
    Kx   = Vector{Tf}(undef,nnzK)

    #fill in triu(P)
    @inbounds for cidx = 1:n
        for j = (P.colptr[cidx]):(P.colptr[cidx+1]-1)
            ridx = P.rowval[j]
            if(ridx <= cidx)
                Kidx = KNextInCol[cidx]
                Ki[Kidx] = ridx
                Kx[Kidx] = P.nzval[j]
                KNextInCol[cidx] += 1
            else
                break
            end
        end
    end

    #fill in A' in the upper RHS
    @inbounds for cidx = 1:n
        for j = (A.colptr[cidx]):(A.colptr[cidx+1]-1)
            ridx = A.rowval[j]
            ridxpn = ridx + n
            Kidx = KNextInCol[ridxpn]
            KNextInCol[ridxpn] += 1
            Ki[Kidx] = cidx
            Kx[Kidx] = A.nzval[j]
        end
    end

    #fill in the Sigma part
    @inbounds for cidx = 1:n
        if(KNextInCol[cidx] == Kp[cidx+1])        #P had diagonal term at cidx)
            Kx[KNextInCol[cidx] - 1] += S[cidx]
        else
            Kx[KNextInCol[cidx]]      = S[cidx]
            Ki[KNextInCol[cidx]]      = cidx
        end
    end
    #fill in the rho part
    @inbounds for idx = 1:m
        kidx = KNextInCol[idx+n]
        Kx[kidx] = R[idx]
        Ki[kidx] = idx+n
    end

    #assemble the matrix
    K = SparseMatrixCSC(m+n,m+n,Kp,Ki,Kx)

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
        K   = _kktutils_make_kkt(P, A, sigma, rho, :U)
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
