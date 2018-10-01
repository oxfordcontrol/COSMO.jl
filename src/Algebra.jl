using LinearAlgebra
import LinearAlgebra; lmul!, rmul!
export colNorms!, rowNorms!, lrmul!
const IdentityMatrix = UniformScaling{Bool}

function colNorms!(v::Array{Tf,1},
    A::Matrix{Tf};
    reset::Bool = true) where{Tf<:AbstractFloat}

    if(reset) v.= 0 end

    for i=1:size(A,2)
        v[i] = max(v[i],norm(view(A,:,i),Inf))
    end
end

function colNorms!(v::Array{Tf,1},
    A::SparseMatrixCSC{Tf,Ti}; reset::Bool = true) where{Tf<:AbstractFloat,Ti<:Integer}

    if(reset) v.= 0 end

    for i=1:A.n
        for j = A.colptr[i]:(A.colptr[i+1] - 1)
            v[i] = max(v[i],abs(A.nzval[j]))
        end
    end
end

function rowNorms!(v::Array{Tf,1},
    A::Matrix{Tf};
    reset::Bool = true) where{Tf<:AbstractFloat}

    if(reset) v.= 0 end

    for i=1:size(A,1)
        v[i] = max(v[i],norm(view(A,i,:),Inf))
    end
end



function rowNorms!(v::Array{Tf,1},
    A::SparseMatrixCSC{Tf,Ti};
    reset::Bool = true) where{Tf<:AbstractFloat,Ti<:Integer}

    if(reset) v.= 0 end

    for i=1:(A.colptr[end]-1)
        v[A.rowval[i]] = max(v[A.rowval[i]],abs(A.nzval[i]))
    end

end


function LinearAlgebra.lmul!(L::IdentityMatrix, M)
    return M
end

function LinearAlgebra.rmul!(M,R::IdentityMatrix)
    return M
end

function LinearAlgebra.lmul!(L::Diagonal, M::SparseMatrixCSC)

    #NB : Same as:  @views M.nzval .*= D.diag[M.rowval]
    #but this way allocates no memory at all and
    #is marginally faster
    for i = 1:(M.colptr[end]-1)
        M.nzval[i] *= L.diag[M.rowval[i]]
    end
    M
end

function LinearAlgebra.rmul!(M::SparseMatrixCSC,R::Diagonal)
    for i = 1:M.n
        for j = M.colptr[i]:(M.colptr[i+1]-1)
            M.nzval[j] *= R.diag[i]
        end
    end
    M
end

function lrmul!(L::IdentityMatrix, M::AbstractMatrix, R::IdentityMatrix)
    return M
end

function lrmul!(L::Diagonal, M::AbstractMatrix, R::Diagonal)
    lmul!(L,rmul!(M,R))
end

function lrmul!(L::Diagonal, M::AbstractMatrix, R::IdentityMatrix)
    lmul!(L,M)
end

function lrmul!(L::IdentityMatrix, M::AbstractMatrix, R::Diagonal)
    rmul!(M,R)
end

function lrmul!(L::Diagonal, M::SparseMatrixCSC, R::Diagonal)
    for i = 1:M.n
        for j = M.colptr[i]:(M.colptr[i+1]-1)
            M.nzval[j] *= L.diag[M.rowval[j]]*R.diag[i]
        end
    end
    M
end
