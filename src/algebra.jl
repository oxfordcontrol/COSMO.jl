using  LinearAlgebra
const  IdentityMatrix = UniformScaling{Bool}


function clip(s::Real, min_thresh::Real, max_thresh::Real, min_new::Real = min_thresh, max_new::Real = max_thresh)
	s = ifelse(s < min_thresh, min_new, ifelse(s > max_thresh, max_new, s))
end

function scaled_norm(E::IdentityMatrix, v::Array, p::Real = 2)
	E.λ ? norm(v, p) : zero(eltype(v))
end

function scaled_norm(E::Diagonal, v::Array{T}, p::Real = 2) where{T}
	if p == 2
		return scaled_norm2(E, v)
	elseif p == Inf
		return scaled_norm_Inf(E, v)
	elseif p == 1
		return scaled_norm1(E, v)
	else
		throw(ArgumentError("bad norm specified"))
	end
end

function scaled_norm2(E::Diagonal, v::Array)
	sum_sq  = zero(eltype(v))
	for i = 1:length(v)
		sum_sq += (E.diag[i] * v[i])^2
	end
	return sqrt(sum_sq)
end

function scaled_norm_Inf(E::Diagonal, v::Array)
	norm  = zero(eltype(v))
	for i = 1:length(v)
		norm = max(norm, abs(E.diag[i] * v[i]))
	end
	return norm
end

function scaled_norm1(E::Diagonal, v::Array)
	norm  = zero(eltype(v))
	for i = 1:length(v)
		norm += abs(E.diag[i] * v[i])
	end
	return norm
end

function col_norms!(v::Array{Tf, 1},
	A::Matrix{Tf};
	reset::Bool = true) where{Tf <: AbstractFloat}

	if(reset) v.= 0 end

	for i = 1:size(A, 2)
		v[i] = max(v[i], norm(view(A, :, i), Inf))
	end
	return v
end

function col_norms!(v::Array{Tf, 1},
	A::SparseMatrixCSC{Tf,Ti}; reset::Bool = true) where{Tf <: AbstractFloat, Ti <: Integer}

	if(reset) v.= 0 end

	for i = 1:A.n
		@inbounds for j = A.colptr[i]:(A.colptr[i + 1] - 1)
		@inbounds v[i] = max(v[i], abs(A.nzval[j]))
	end
	end
	return v
end

function row_norms!(v::Array{Tf, 1},
	A::Matrix{Tf};
	reset::Bool = true) where{Tf <: AbstractFloat}

	if(reset) v.= 0 end

	for i = 1:size(A, 1)
		v[i] = max(v[i], norm(view(A, i, :), Inf))
	end
	return v
end

function row_norms!(v::Array{Tf, 1},
	A::SparseMatrixCSC{Tf, Ti};
	reset::Bool = true) where{Tf <: AbstractFloat, Ti <: Integer}

	if(reset) v.= 0 end

	@inbounds for i = 1:(A.colptr[end] - 1)
		@inbounds v[A.rowval[i]] = max(v[A.rowval[i]], abs(A.nzval[i]))
	end
	return v
end

function lmul!(L::Diagonal, M::SparseMatrixCSC)

	#NB : Same as:  @views M.nzval .*= D.diag[M.rowval]
	#but this way allocates no memory at all and
	#is marginally faster
	m, n = size(M)
	(m == length(L.diag)) || throw(DimensionMismatch())

	@inbounds for i = 1:(M.colptr[end] - 1)
	@inbounds M.nzval[i] *= L.diag[M.rowval[i]]
	end
	return M
end

lmul!(L::IdentityMatrix, M::AbstractMatrix) = L.λ ? M : M .= zero(eltype(M))

function rmul!(M::SparseMatrixCSC, R::Diagonal)

	m, n = size(M)
	(n == length(R.diag)) || throw(DimensionMismatch())

	@inbounds for i = 1:n, j = M.colptr[i]:(M.colptr[i + 1] - 1)
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
	(m == length(Ld) && n == length(Rd)) || throw(DimensionMismatch())

	@inbounds for i = 1:n, j = Mcolptr[i]:(Mcolptr[i + 1] - 1)
	@inbounds Mnzval[j] *= Ld[Mrowval[j]] * Rd[i]
	end
	return M
end

lrmul!(L::IdentityMatrix,
	M::AbstractMatrix,
	R::IdentityMatrix) = (L.λ && R.λ) ? M : M .= zero(eltype(M))

lrmul!(L::Diagonal,
	M::SparseMatrixCSC,
	R::IdentityMatrix) = R.λ ? lmul!(L, M) : M .= zero(eltype(M))

lrmul!(L::Diagonal,
	M::AbstractMatrix,
	R::Diagonal) = LinearAlgebra.lmul!(L, LinearAlgebra.rmul!(M, R))

lrmul!(L::Diagonal,
	M::AbstractMatrix,
	R::IdentityMatrix) = R.λ ? LinearAlgebra.lmul!(L, M) : M .= zero(eltype(M))


lrmul!(L::IdentityMatrix,
	M::AbstractMatrix,
	R::Diagonal) = L.λ ? LinearAlgebra.rmul!(M, R) : M .= zero(eltype(M))

"""
    symmetrize_upper!(A)

Symmetrizes the matrix A by calculating A = 0.5 * (A + A') but only performs the operation on the upper triangular part.
"""
function symmetrize_upper!(A::AbstractMatrix)
	n = size(A, 1)
	@assert(size(A, 1) == size(A, 2), "Matrix is not square.")
	@inbounds for j in 1:n, i in 1:j
 		A[i, j] = (A[i, j] + A[j, i]) / 2
	end
	nothing
end

"""
    symmetrize_full!(A)

Symmetrizes the matrix A by calculating A = 0.5 * (A + A') and storing the result in-place.
"""
function symmetrize_full!(A::AbstractMatrix)
	n = size(A, 1)
	@assert(size(A, 1) == size(A, 2), "Matrix is not square.")
	@inbounds for j in 1:n, i in 1:j
 		A[i, j] = (A[i, j] + A[j, i]) / 2
 		A[j, i] = A[i, j]
	end
	nothing
end

# this function assumes real symmetric X and only considers the upper triangular part
function is_pos_sem_def(X, tol)
    # set option 'N' to only compute eigenvalues, s is ordered from min to max
   s, U = LAPACK.syevr!('N', 'A', 'U', X, 0.0, 0.0, 0, 0, -1.0);
   return s[1] >= -tol
end

function is_neg_sem_def(X, tol)
    # set option 'N' to only compute eigenvalues, s is ordered from min to max
   s, U = LAPACK.syevr!('N', 'A', 'U', X, 0.0, 0.0, 0, 0, -1.0);
   return s[end] <= tol
end
