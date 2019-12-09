using  LinearAlgebra
import LinearAlgebra.lmul!
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

function lmul!(L::Diagonal, x::AbstractVector)	
	(length(L.diag) == length(x)) || throw(DimensionMismatch())	
	@. x = x * L.diag	
	return nothing	
end	

lmul!(L::IdentityMatrix, x::AbstractVector) = L.λ ? x : x .= zero(eltype(x))

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
function is_pos_def!(X::AbstractMatrix{T}, tol::T=zero(T)) where T
	# See https://math.stackexchange.com/a/13311
	@inbounds for i = 1:size(X, 1)
		X[i, i] += tol
	end
	F = cholesky!(Symmetric(X), check = false)
	return issuccess(F)
end

function is_neg_def!(X::AbstractMatrix{T}, tol::T=zero(T)) where T
	@. X *= -one(T)
	return is_pos_def!(X, tol)
end

is_pos_def(X::AbstractMatrix{T}, tol::T=zero(T)) where T = is_pos_def!(copy(X), tol)
is_neg_def(X::AbstractMatrix{T}, tol::T=zero(T)) where T = is_pos_def!(-X, tol)

function populate_upper_triangle!(A::AbstractMatrix{T}, x::AbstractVector{T}, scaling_factor::T=T(1/sqrt(2))) where T
	k = 0
	 for j in 1:size(A, 2)
		for i in 1:j
		   k += 1
		   if i != j
			   A[i, j] = scaling_factor * x[k]
		   else
			   A[i, j] = x[k]
		   end
		 end
	 end
	 nothing
end

function extract_upper_triangle!(A::AbstractMatrix{T}, x::AbstractVector{T}, scaling_factor::T=T(sqrt(2))) where T
   k = 0
	 for j in 1:size(A, 2)
		for i in 1:j
		   k += 1
		   if i != j
			   x[k] = scaling_factor * A[i, j]
		   else
			   x[k] = A[i, j]
		   end
		 end
	 end
   nothing
end

function populate_upper_triangle(x::AbstractVector{T}, scaling_factor::T=T(1/sqrt(2))) where T
	n = Int(1/2*(sqrt(8*length(x) + 1) - 1)) # Solution of (n^2 + n)/2 = length(x) obtained by WolframAlpha
	A = zeros(n, n)
	populate_upper_triangle!(A, x, scaling_factor)
	return Symmetric(A)
end

function extract_upper_triangle(A::AbstractMatrix{T}, scaling_factor::T=T(sqrt(2))) where T
	n = size(A, 2)
	x = zeros(Int(n*(n + 1)/2))
	extract_upper_triangle!(A, x, scaling_factor)
	return x
end