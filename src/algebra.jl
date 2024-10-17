using  LinearAlgebra
const  IdentityMatrix = UniformScaling{Bool}


function clip(s::Real, min_thresh::Real, max_thresh::Real, min_new::Real = min_thresh, max_new::Real = max_thresh)
	s = ifelse(s < min_thresh, min_new, ifelse(s > max_thresh, max_new, s))
end

function scaled_norm(E::IdentityMatrix, v::Array, p::T = T(2)) where {T <: AbstractFloat}
	E.λ ? norm(v, p) : zero(eltype(v))
end

function scaled_norm(E::Diagonal{T}, v::Array{T}, p::T = T(2)) where {T <: AbstractFloat}
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
	reset::Bool = true) where {Tf <: AbstractFloat}

	if reset
		fill!(v,0.)
	end

	for i = 1:size(A, 2)
		v[i] = max(v[i], norm(view(A, :, i), Inf))
	end
	return v
end

function col_norms!(v::Array{Tf, 1},
	A::SparseMatrixCSC{Tf,Ti}; reset::Bool = true) where {Tf <: AbstractFloat, Ti <: Integer}

	if reset
		fill!(v, 0)
	end

	@inbounds for i = eachindex(v)
		for j = A.colptr[i]:(A.colptr[i + 1] - 1)
			tmp = abs(A.nzval[j])
			v[i] = v[i] > tmp ? v[i] : tmp;
		end
	end
	return v
end

function row_norms!(v::Array{Tf, 1},
	A::Matrix{Tf};
	reset::Bool = true) where{Tf <: AbstractFloat}

	if reset
		fill!(v,0.)
	end

	for i = 1:size(A, 1)
		v[i] = max(v[i], norm(view(A, i, :), Inf))
	end
	return v
end

function row_norms!(v::Array{Tf, 1},
	A::SparseMatrixCSC{Tf, Ti};
	reset::Bool = true) where{Tf <: AbstractFloat, Ti <: Integer}

	if reset
		fill!(v,0.)
	end

	@inbounds for i = 1:(A.colptr[end] - 1)
		idx = A.rowval[i]
		tmp = abs(A.nzval[i])
		v[idx] = v[idx] > tmp ? v[idx] : tmp
	end
	return v
end

function scalarmul!(A::SparseMatrixCSC, c::Real)
	A.nzval .*= c
end

function scalarmul!(A::AbstractMatrix, c::Real)
	A .*= c
end


function lmul!(L::Diagonal{T}, M::SparseMatrixCSC{T}) where {T <: AbstractFloat}

	#NB : Same as:  @views M.nzval .*= D.diag[M.rowval]
	#but this way allocates no memory at all and
	#is marginally faster
	m, n = size(M)
	(m == length(L.diag)) || throw(DimensionMismatch())

	@inbounds for i = 1:(M.colptr[end] - 1)
	 		M.nzval[i] *= L.diag[M.rowval[i]]
	end
	return M
end

lmul!(L::IdentityMatrix, M::AbstractMatrix) = L.λ ? M : M .= zero(eltype(M))

function lmul!(L::Diagonal{T}, x::AbstractVector{T}) where {T <: AbstractFloat}
	(length(L.diag) == length(x)) || throw(DimensionMismatch())
	@. x = x * L.diag
	return nothing
end

lmul!(L::IdentityMatrix, x::AbstractVector{T}) where {T <: AbstractFloat} = L.λ ? x : x .= zero(eltype(x))



function rmul!(M::SparseMatrixCSC{T}, R::Diagonal{T}) where {T <: AbstractFloat}

	m, n = size(M)
	(n == length(R.diag)) || throw(DimensionMismatch())

	@inbounds for i = 1:n, j = M.colptr[i]:(M.colptr[i + 1] - 1)
		 	M.nzval[j] *= R.diag[i]
	end
	return M
end

rmul!(M::AbstractMatrix, R::IdentityMatrix) = R.λ ? R : R .= zero(eltype(R))

function lrmul!(L::Diagonal{T}, M::SparseMatrixCSC{T}, R::Diagonal{T}) where {T <: AbstractFloat}

	m, n = size(M)
	Mnzval  = M.nzval
	Mrowval = M.rowval
	Mcolptr = M.colptr
	Rd      = R.diag
	Ld      = L.diag
	(m == length(Ld) && n == length(Rd)) || throw(DimensionMismatch())

	@inbounds for i = 1:n
		for j = Mcolptr[i]:(Mcolptr[i + 1] - 1)
	 		Mnzval[j] *= Ld[Mrowval[j]] * Rd[i]
		end
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
function symmetrize_upper!(A::AbstractMatrix{T}) where {T <: AbstractFloat}
	n = size(A, 1)
	@assert(size(A, 1) == size(A, 2), "Matrix is not square.")
	@inbounds for j in 1:n, i in 1:j
 		A[i, j] = (A[i, j] + A[j, i]) / T(2)
	end
	nothing
end

"""
    symmetrize_full!(A)

Symmetrizes the matrix A by calculating A = 0.5 * (A + A') and storing the result in-place.
"""
function symmetrize_full!(A::AbstractMatrix{T}) where {T <: AbstractFloat}
	n = size(A, 1)
	@assert(size(A, 1) == size(A, 2), "Matrix is not square.")
	@inbounds for j in 1:n, i in 1:j
 		A[i, j] = (A[i, j] + A[j, i]) / T(2)
 		A[j, i] = A[i, j]
	end
	nothing
end

# this function assumes real or complex Hermitian X and only considers the upper triangular part
function is_pos_def!(X::AbstractMatrix{<:RealOrComplex{T}}, tol::T=zero(T)) where {T}
	# See https://math.stackexchange.com/a/13311
	@inbounds for i = 1:size(X, 1)
		X[i, i] += tol
	end
	F = cholesky!(Hermitian(X), check = false)
	return issuccess(F)
end

function is_neg_def!(X::AbstractMatrix{<:RealOrComplex{T}}, tol::T=zero(T)) where {T}
	@. X *= -one(T)
	return is_pos_def!(X, tol)
end

is_pos_def(X::AbstractMatrix{<:RealOrComplex{T}}, tol::T=zero(T)) where {T} = is_pos_def!(copy(X), tol)
is_neg_def(X::AbstractMatrix{<:RealOrComplex{T}}, tol::T=zero(T)) where {T} = is_pos_def!(-X, tol)


"Round x to the closest multiple of N."
function round_multiple(x::T, N::T) where {T <: Integer}
	return floor(T, x +  0.5 * N - rem(x + 0.5 * N, N))
end
