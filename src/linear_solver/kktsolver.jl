export QdldlKKTSolver, CholmodKKTSolver
# -------------------------------------
# abstract type defs
# -------------------------------------
abstract type AbstractKKTSolver end
abstract type AbstractKKTShape end

# NB: all concrete examples of this type should
# implement refactor! and solve! methods and should
# implement a constructor taking (P,A,sigma,rho) as
# arguments

# -------------------------------------
# some internal utility functions
# -------------------------------------

"Check the dimensions of the problem data used to assemble the KKT matrix `K`."
function _kktutils_check_dims(P::SparseMatrixCSC{Tv, Ti}, A::AbstractMatrix{Tv}, sigma::Tv, rho::Union{Tv, AbstractVector{Tv}}) where {Tv <: AbstractFloat, Ti <: Integer}

    n = size(P, 1)
    m = size(A, 1)

    size(A,2) == n || throw(DimensionMismatch())

    length(rho)   == m || length(rho)   == 1 || throw(DimensionMismatch())
    length(sigma) == n || length(sigma) == 1 || throw(DimensionMismatch())

    return m, n
end

"Count the number of nnz in each column of the lower triangle of the KKT matrix `K`."
function _count_lower_triangle!(Kcolnz::Vector{Ti}, P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, n::Int) where {Tv <: AbstractFloat, Ti <: Integer}
    # count strict lower triangular values of P
    @inbounds for cidx = 1:n
        for j = (P.colptr[cidx]):(P.colptr[cidx+1]-1)
            ridx = P.rowval[j]
            if (ridx > cidx)
                Kcolnz[cidx] += 1
            end
        end
    end
    # add number of column entries of A
    @inbounds for cidx = 1:n
        Kcolnz[cidx] += A.colptr[cidx+1] - A.colptr[cidx]
    end

    #every element on the diagonal is also nonzero
    Kcolnz .+= 1
end

"""
Count the number of nnz in each column of the upper triangle of the KKT matrix `K`.
"""
function _count_upper_triangle!(Kcolnz::Vector{<:Integer}, P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, n::Integer) where {Tv <: AbstractFloat, Ti <: Integer}
    # count nonzeros in the strict upper triangle of (P)
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
end

"Determine and fill in the rowvals `Ki` and nzvals `Kx` of the lower triangle of `K`."
function _fill_lower_triangle!(Ki::Vector{Ti}, Kp::Vector{Ti}, Kx::Vector{Tv}, P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, R::Union{Tv, AbstractVector{Tv}}, S::Union{Tv, AbstractVector{Tv}}, KNextInCol::Vector{Ti}, m::Int, n::Int) where {Tv <: AbstractFloat, Ti <: Integer}
    # tril(K) = | tril(P) + σI |         |
    #           |   A          | - 1/ρ I |


    #fill in the Sigma part
    @inbounds for cidx = 1:n
            Kx[KNextInCol[cidx]]      = S[cidx]
            Ki[KNextInCol[cidx]]      = cidx
            KNextInCol[cidx] += 1
    end

    #fill in tril(P)
    @inbounds for cidx = 1:n
        for j = (P.colptr[cidx]):(P.colptr[cidx+1]-1)
            ridx = P.rowval[j]
            if ridx > cidx
                Kidx = KNextInCol[cidx]
                Ki[Kidx] = ridx
                Kx[Kidx] = P.nzval[j]
                KNextInCol[cidx] += 1
            elseif ridx == cidx # add diagonal entries of P onto the sigma value
                Kidx = KNextInCol[cidx] - 1
                Kx[Kidx] += P.nzval[j]
            end
        end
    end

    #fill in A in the lower LHC
    @inbounds for cidx = 1:n
        for j = (A.colptr[cidx]):(A.colptr[cidx+1]-1)
            ridx = A.rowval[j]
            ridxpn = ridx + n
            Kidx = KNextInCol[cidx]
            Ki[Kidx] = ridxpn
            Kx[Kidx] = A.nzval[j]
            KNextInCol[cidx] += 1
        end
    end

    #fill in the rho part
    @inbounds for idx = 1:m
        kidx = KNextInCol[idx + n]
        Kx[kidx] = R[idx]
        Ki[kidx] = idx + n
    end
end

"Determine and fill in the rowvals `Ki` and nzvals `Kx` of the upper triangle of `K`."
function _fill_upper_triangle!(Ki::Vector{Ti}, Kp::Vector{Ti}, Kx::Vector{Tv}, P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, R::Union{Tv, AbstractVector{Tv}}, S::Union{Tv, AbstractVector{Tv}}, KNextInCol::Vector{Ti}, m::Int, n::Int) where {Tv <: AbstractFloat, Ti <: Integer}
    # triu(K) = | triu(P) + σI |    A'   |
    #           |              | - 1/ρ I |

    #fill in triu(P)
    @inbounds for cidx = 1:n
        for j = (P.colptr[cidx]):(P.colptr[cidx+1]-1)
            ridx = P.rowval[j]
            if ridx <= cidx
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
        if KNextInCol[cidx] == Kp[cidx+1]       #P had diagonal term at cidx)
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
end

# handle the different forms of inputs for sigma and rho: here both are vectors (*)
function assemble_kkt_triangle(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma::AbstractVector{Tv}, rho::AbstractVector{Tv}, shape::Symbol = :U) where {Tv <: AbstractFloat, Ti <: Integer}
        S = sigma
        R  = Array{Tv}(undef, size(A, 1))
     @. R  = -one(Tv) / rho
     # pass on to main function
    return _assemble_kkt_triangle(P, A, S, R, shape)
end
# only sigma is a scalar
function assemble_kkt_triangle(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma::Tv, rho::AbstractVector{Tv}, shape::Symbol = :U) where {Tv <: AbstractFloat, Ti <: Integer}
    S = Vector{Tv}(undef, size(P, 1))
    fill!(S, sigma)
    # put again through (*)
    return assemble_kkt_triangle(P, A, S, rho, shape)
end

# only rho is a scalar
function assemble_kkt_triangle(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma::AbstractVector{Tv}, rho::Tv, shape::Symbol = :U) where {Tv <: AbstractFloat, Ti <: Integer}
    R  = Array{Tv}(undef, size(A, 1))
    fill!(R, -one(Tv) ./ rho)
    # pass on to main function
    return _assemble_kkt_triangle(P, A, sigma, R, shape)
end

# rho and sigma are scalars
function assemble_kkt_triangle(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma::Tv, rho::Tv, shape::Symbol = :U) where {Tv <: AbstractFloat, Ti <: Integer}
    S = Vector{Tv}(undef, size(P, 1))
    fill!(S, sigma)

    R  = Array{Tv}(undef, size(A, 1))
    fill!(R, -one(Tv) ./ rho)
    # pass on to main function
    return _assemble_kkt_triangle(P, A, S, R, shape)
end
"Given `P`, `A`, `sigma` and `rho` return the upper / lower triangle, defined by `shape`,  of the KKT condition matrix `K`."
function _assemble_kkt_triangle(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, S::AbstractVector{Tv}, R::AbstractVector{Tv}, shape::Symbol = :U) where {Tv <: AbstractFloat, Ti <: Integer}

    #   K = | P + σI |  A'     |
    #       |   A    | - 1/ρ I |
    #         left      right
    n = size(P, 1)
    m = size(A, 1)

    # make the column counter of K
    Kcolnz = zeros(Ti, m + n)

    #work out how many nonzeros are in K.
    if shape == :U
        _count_upper_triangle!(Kcolnz, P, A, n)
    elseif shape == :L
        _count_lower_triangle!(Kcolnz, P, A, n)
    end

    #make the column pointer and nonzero count
    Kp   = Vector{Ti}(undef, m + n + 1)
    Kp[1] = 1
    @inbounds for i = 1:(m+n)
        Kp[i+1] = Kp[i] + Kcolnz[i]
    end

    KNextInCol = copy(Kp);        #next value in column goes here
    nnzK = Kp[end] - 1
    Ki   = Vector{Ti}(undef, nnzK)
    Kx   = Vector{Tv}(undef, nnzK)

    # fill in the rowvals and nzvals
    if shape == :U
        _fill_upper_triangle!(Ki, Kp, Kx, P, A, R, S, KNextInCol, m, n)
    elseif shape == :L
        _fill_lower_triangle!(Ki, Kp, Kx, P, A, R, S, KNextInCol, m, n)
    end

    #assemble the matrix
    K = SparseMatrixCSC(m + n, m + n, Kp, Ki, Kx)

    return K
end

"Given `P`, `A`, `sigma` and `rho` return the full KKT condition matrix `K`."
function _assemble_kkt_full(P::SparseMatrixCSC{Tv, Ti}, A::AbstractMatrix{Tv}, sigma::Tv, rho::Union{Tv, AbstractVector{Tv}}) where {Tv <: AbstractFloat, Ti <: Integer}

    n = size(P, 1)
    m = size(A, 1)
    S = length(sigma) == 1 ? (sigma[1]) * I : Diagonal(sigma)
    rhoinv  = Vector{Tv}(undef, m)
    rhoinv .= (-one(Tv)./rho)
    D       = SparseMatrixCSC(Diagonal(rhoinv))

    #compute the full KKT matrix
    K = [P + S SparseMatrixCSC(A'); A D]

    return K
end

"Update the diagonal `rho` entries in the lower right hand corner of the KKT matrix `K`."
function update_kkt_matrix!(K::SparseMatrixCSC{Tv, Ti}, n::Int, m::Int, rho::Tv) where {Tv <: AbstractFloat, Ti <: Integer}
    @inbounds @simd for i = (n + 1):(n + m)
        K[i, i] = -one(Tv) / rho
    end
end

function update_kkt_matrix!(K::SparseMatrixCSC{Tv, Ti}, n::Int, m::Int, rho_vec::AbstractVector{Tv}) where {Tv <: AbstractFloat, Ti <: Integer}
    @inbounds @simd for i = (n + 1):(n + m)
        K[i, i] = -one(Tv) / rho_vec[i - n]
    end
end

# -------------------------------------
# QDLDL solver
# -------------------------------------

struct QdldlKKTSolver{Tv, Ti} <: AbstractKKTSolver

    m::Ti
    n::Ti
    ldlfact::QDLDL.QDLDLFactorisation{Tv, Ti}
    diagidx::Vector{Ti}

    function QdldlKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma::Tv, rho::Union{Tv, AbstractVector{Tv}}) where {Tv <: AbstractFloat, Ti <: Integer}
        m, n = _kktutils_check_dims(P, A, sigma, rho)
        K    = assemble_kkt_triangle(P, A, sigma, rho, :U)
        ldlfact = qdldl(K)
        
        #we want the index if diagonal entries with K triu, 
        #this means we want the *last* entry in every column
        #and assume the diagonal is fully populated
        @views diagidx = K.colptr[2:end] .- 1
        

        #check for exactly n positive eigenvalues
        positive_inertia(ldlfact) == n || error("Objective function is not convex.")

        new{Tv, Ti}(m, n, ldlfact,diagidx)
    end
end

function solve!(s::QdldlKKTSolver, lhs, rhs)
    @. lhs = rhs
    QDLDL.solve!(s.ldlfact, lhs)
end


function update_rho!(s::QdldlKKTSolver{Tv, Ti}, rho::Union{Tv, AbstractVector{Tv}}) where {Tv <: AbstractFloat, Ti <: Integer}

    @views QDLDL.update_values!(s.ldlfact, s.diagidx[(s.n+1):(s.n+s.m)], (-one(Ti) ./ rho))
    refactor!(s.ldlfact)
end

# -------------------------------------
# Julia Native solver (CHOLMOD based)
# -------------------------------------
mutable struct CholmodKKTSolver{Tv, Ti} <: AbstractKKTSolver

    fact::SuiteSparse.CHOLMOD.Factor{Float64} # SuiteSparse's cholmod only supports doubles
    K::SparseMatrixCSC{Tv, Ti}
    m::Ti
    n::Ti

    function CholmodKKTSolver(P::SparseMatrixCSC{Tv, Ti}, A::SparseMatrixCSC{Tv, Ti}, sigma::Tv, rho) where {Tv <: AbstractFloat, Ti <: Integer}

         m,n  = _kktutils_check_dims(P, A, sigma, rho)
         K    = _assemble_kkt_full(P, A, sigma, rho)
         fact = ldlt(K)

        return new{Tv, Ti}(fact, K, m, n)

    end
end

solve!(s::CholmodKKTSolver, lhs, rhs) = lhs .= s.fact \ rhs

function update_rho!(s::CholmodKKTSolver{Tv, Ti}, rho::Union{Tv, AbstractVector{Tv}}) where {Tv <: AbstractFloat, Ti <: Integer}
    update_kkt_matrix!(s.K, s.n, s.m, rho)
    #complete restart
    s.fact = ldlt(s.K)
end

free_memory!(s::AbstractKKTSolver) = nothing
