using  LinearAlgebra
import LinearAlgebra: lmul!, rmul!
const  IdentityMatrix = UniformScaling{Bool}


function clip(s,minThresh,maxThresh,minNew = minThresh,maxNew = maxThresh)
    s = ifelse(s < minThresh,minNew,ifelse(s > maxThresh,maxNew,s))
end

function scalednorm(E::IdentityMatrix,v::Array,p::Real=2)
    E.λ ? norm(v,p) : zero(eltype(v))
end

function scalednorm(E::Diagonal,v::Array{T},p::Real=2) where{T}
    if p == 2
        return scalednorm2(E,v)::T
    elseif p == Inf
        return scalednormInf(E,v)::T
    elseif p == 1
        return scaledNorm1(E,v)::T
    else
        throw(ArgumentError("bad norm specified"))
    end
end

function scalednorm2(E::Diagonal,v::Array)

    sumsq  = zero(eltype(v))
    for i = 1:length(v)
        sumsq += (E.diag[i]*v[i])^2
    end
    return sqrt(sumsq)::eltype(v)
end

function scalednormInf(E::Diagonal,v::Array)
    norm  = zero(eltype(v))
    for i = 1:length(v)
        norm = max(norm,(E.diag[i]*v[i]))
    end
    return norm::eltype(v)
end

function scalednorm1(E::Diagonal,v::Array)
    norm  = zero(eltype(v))
    for i = 1:length(v)
        norm += abs(E.diag[i]*v[i])
    end
    return norm::eltype(v)
end



function colnorms!(v::Array{Tf,1},
    A::Matrix{Tf};
    reset::Bool = true) where{Tf<:AbstractFloat}

    if(reset) v.= 0 end

    for i=1:size(A,2)
        v[i] = max(v[i],norm(view(A,:,i),Inf))
    end
    return v
end

function colnorms!(v::Array{Tf,1},
    A::SparseMatrixCSC{Tf,Ti}; reset::Bool = true) where{Tf<:AbstractFloat,Ti<:Integer}

    if(reset) v.= 0 end

    for i=1:A.n
        @inbounds for j = A.colptr[i]:(A.colptr[i+1] - 1)
            @inbounds v[i] = max(v[i],abs(A.nzval[j]))
        end
    end
    return v
end

function rownorms!(v::Array{Tf,1},
    A::Matrix{Tf};
    reset::Bool = true) where{Tf<:AbstractFloat}

    if(reset) v.= 0 end

    for i=1:size(A,1)
        v[i] = max(v[i],norm(view(A,i,:),Inf))
    end
    return v
end

function rownorms!(v::Array{Tf,1},
    A::SparseMatrixCSC{Tf,Ti};
    reset::Bool = true) where{Tf<:AbstractFloat,Ti<:Integer}

    if(reset) v.= 0 end

    @inbounds for i=1:(A.colptr[end]-1)
        @inbounds v[A.rowval[i]] = max(v[A.rowval[i]],abs(A.nzval[i]))
    end
    return v
end

function lmul!(L::Diagonal, M::SparseMatrixCSC)

    #NB : Same as:  @views M.nzval .*= D.diag[M.rowval]
    #but this way allocates no memory at all and
    #is marginally faster
    m, n = size(M)
    (m==length(L.diag)) || throw(DimensionMismatch())

    @inbounds for i = 1:(M.colptr[end]-1)
        @inbounds M.nzval[i] *= L.diag[M.rowval[i]]
    end
    return M
end

lmul!(L::IdentityMatrix, M::AbstractMatrix) = R.λ ? R : R .= zero(eltype(M))

function rmul!(M::SparseMatrixCSC,R::Diagonal)

    m, n = size(M)
    (n==length(R.diag)) || throw(DimensionMismatch())

    @inbounds for i = 1:n, j = M.colptr[i]:(M.colptr[i+1]-1)
        @inbounds M.nzval[j] *= R.diag[i]
    end
    return M
end

rmul!(M::AbstractMatrix, R::IdentityMatrix) = R.λ ? R : R .= zero(eltype(R))

function lrmul!(L::Diagonal, M::SparseMatrixCSC, R::Diagonal)

    Mnzval  = M.nzval
    Mrowval = M.rowval
    Mcolptr = M.colptr
    Rd      = R.diag
    Ld      = L.diag

    m, n = size(M)
    (m==length(Ld) && n==length(Rd)) || throw(DimensionMismatch())

    @inbounds for i = 1:n, j = Mcolptr[i]:(Mcolptr[i+1]-1)
        @inbounds Mnzval[j] *= Ld[Mrowval[j]]*Rd[i]
    end
    return M
end

lrmul!(L::IdentityMatrix,
       M::AbstractMatrix,
       R::IdentityMatrix) = (L.λ && R.λ) ? M : M .= zero(eltype(M))

lrmul!(L::Diagonal,
       M::AbstractMatrix,
       R::Diagonal) = lmul!(L,rmul!(M,R))

lrmul!(L::Diagonal,
       M::AbstractMatrix,
       R::IdentityMatrix) = R.λ ? lmul!(L,M) : M .= zero(eltype(M))

lrmul!(L::IdentityMatrix,
       M::AbstractMatrix,
       R::Diagonal      ) = L.λ ? rmul!(M,R) : M .= zero(eltype(M))
