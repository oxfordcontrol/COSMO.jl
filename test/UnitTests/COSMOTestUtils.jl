
# generate a random pos def matrix with eigenvalues between 0.1 and 2
function generate_pos_def_matrix(rng::MersenneTwister, n::Int64, aMin::Real = 0.1, aMax::Real = 2)
	X = rand(rng, n, n)
	# any real square matrix can be QP decomposed into a orthogonal matrix and an uppertriangular matrix R
	Q, R = qr(X)
	eigs = rand(rng ,n) .* (aMax .- aMin) .+ aMin
	X = Q * Matrix(Diagonal(eigs)) * Q'
	X = 0.5 * (X + X')
	return X
end

function is_numerically_pos_sem_def(X, atol)
	X = X ./ 2
	X = X + X'

	F = eigfact(X)
	if size(find( x-> x < -atol, F[:values]), 1) == 0
		return true
	else
		return false
	end
end

function is_numerically_symmetric(X,atol)
	n = size(X, 2)
	for i = 1:n-1, j = i+1:n
		if abs(X[i, j] - X[j, i]) >= atol
			return false
		end
	end
	return true
end


function find_nonsymmetric_component(X)
	for i = 2:size(X, 1), j = 1:(i - 1)
		if abs(X[i, j] - X[j, i]) > 0.0
			return i, j, abs(X[i, j] - X[j, i])
		end
	end
end

function find_different_elements(A, B)
	if size(A) != size(B)
		error("Matrices are not the same size")
	end
	m, n = size(A)
	diff_el = Array[]
	for iii = 1:m, jjj = 1:n
		if A[iii, jjj] != B[iii, jjj]
			push!(diff_el, [iii, jjj])
		end
	end
	return diff_el
end

function duplicate_sparsity_pattern(A)
	m, n = size(A)
	B = zeros(m, n)
	for iii = 1:m, jjj = 1:n
		if A[iii, jjj] != 0
			B[iii, jjj] = 1
		end
	end
	return B
end

function apply_pattern!(A::AbstractMatrix, pattern::AbstractMatrix)
  @assert size(A) == size(pattern) "Matrix A and pattern must have same dimensions."
  m, n = size(A)
  for i = 1:m, j = 1:n
    if pattern[i, j] == 0
      A[i, j] = 0.
    end
  end
end


# create a feasible SDP with one PSDConeTriangle constraint
# min c' x
# s.t. At x + s == bt
#       s ∈ PSDConeTriangle
# choose At and bt in such a way that S has the provided sparsity pattern
function feasible_sdp_with_pattern(rng::MersenneTwister, pattern::AbstractMatrix)
  n = size(pattern, 1)
  d = div(n * (n + 1), 2)


  S = generate_pos_def_matrix(rng, n, 0.1, 2)
  apply_pattern!(S, pattern)
  S = Symmetric(S, :U)

  A1 = rand(rng, n, n)
  apply_pattern!(A1, pattern)
  A1 = Symmetric(A1, :U)
  A = hcat(A1[:])
  s = S[:]
  x = rand(rng, 1)
  b = A * x + s
  B = reshape(b, n, n)

  Y = generate_pos_def_matrix(rng, n,  0.1, 1)
  y = vec(Y)
  P = sparse(zeros(1, 1))
  q = -P * x - A' * y

  Ct = [COSMO.PsdConeTriangle(d)];
  At = zeros(d)
  bt = zeros(d)
  COSMO.extract_upper_triangle!(A1, At, sqrt(2))
  COSMO.extract_upper_triangle!(B, bt, sqrt(2))
  At = hcat(At)
  return P, q, At, bt, Ct
end



"Take a vector `svec` representing the `d` upper-triangular entries of a matrix `X`, and return `X`."
function matrixify(svec::AbstractVector)
  n = - 0.5 + sqrt(0.25 + 2 * length(svec))
  n = Int64(n)
  X = zeros(n, n)
  COSMO.populate_upper_triangle!(X, svec, 1 /sqrt(2))
  return Symmetric(X, :U)
end


function recreate_sparse_matrix(A::SparseMatrixCSC)
	rowInd = A.rowval
	colPtr = A.colptr
	val = A.nzval

	#compute column indices
	colInd = zeros(Int64, length(rowInd))
	cval = 1
	for iii = 2:length(colPtr)
		currentPtr = colPtr[iii]
		prevPtr = colPtr[iii - 1]
		colInd[prevPtr:currentPtr - 1] = cval
		cval += 1
	end

	# sort rowInd and vals
	p = sortperm(rowInd)
	return sparse(rowInd, colInd, val, size(A, 1), size(A, 2))
end

# Geometric mean from https://github.com/JuliaStats/StatsBase.jl
function gmean(a::AbstractArray{T}) where T<:Real
	s = 0.0
	n = length(a)
	for i in 1:n
		tmp = a[i]
		if tmp < 0.0
			throw(DomainError())
		elseif tmp == 0.0
			return 0.0
		else
			s += log(tmp)
		end
	end
	return exp(s / n)
end
