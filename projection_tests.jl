include("src/QOCS.jl")
using Main.QOCS, Test, LinearAlgebra
using BenchmarkTools
using JLD2
using Arpack
using LinearMaps

n = 1000;
X = Symmetric(randn(n, n)) - 40*I; @save "data.jld" X
# @load "data.jld" X 

cone = QOCS.PositiveSemidefiniteCone()
cone.dim = n^2
cone.subspace = zeros(n, n)

function  vecnorm(A) 
	return mapslices(norm, A, dims=1)
end

function my_projection(X, Z, λ, v_rem)
	n = size(X, 1)
	XZ = X*Z
	ZXZ= Z'*XZ

	###
	# Don't propagate parts of the subspace that have already "converged"
	res_norms = vecnorm(XZ - Z*Diagonal(λ))[:]
	sorted_idx = sortperm(res_norms)
	start_idx = findfirst(res_norms[sorted_idx] .> 1e-5)
	if isa(start_idx, Nothing) || start_idx - length(res_norms) > sqrt(n)/2
		start_idx = Int(floor(length(res_norms) - sqrt(n)/2))
	end
	###
	Q = Array(qr(XZ[:, sorted_idx[start_idx:end]] - Z*ZXZ[:, sorted_idx[start_idx:end]]).Q)
	QXZ = Q'*XZ
	XQ = X*Q
	W = Symmetric([ZXZ QXZ'; QXZ Q'*XQ])

	l, V = eigen!(W);
	sorted_idx = sortperm(l)
	idx = findfirst(l[sorted_idx] .> 0) # First positive index
	
	# Positive Ritz pairs
	V_ = V[:, sorted_idx[idx:end]]
	λ = l[sorted_idx[idx:end]]; U = [Z Q]*V_

	# Negative Ritz pairs that we will keep as buffer
	buffer_size = 3;
	Ub = [Z Q]*V[:, sorted_idx[max(idx-buffer_size,1):max(idx-1,1)]]  # Buffer ritz vectors & values
	λb = l[sorted_idx[max(idx-buffer_size,1):max(idx-1,1)]]

	# Projection
	Xπ = U*Diagonal(λ)*U';

	# Estimate largest eigenvalue on the subspace we discarded
	# Careful, we need to enforce all iterates to be orthogonal to the range of U
	# This can be violated by finite precision if we are not careful
	w = zeros(size(U, 2))
	function custom_mul(y::AbstractVector, x::AbstractVector)
		#=
		tmp = Symmetric(X)*x
		y .= tmp - U*(U'*tmp)
		=#
		BLAS.gemv!('T', 1.0, U, x, 0.0, w)
		BLAS.gemv!('N', -1.0, U, w, 1.0, x)
		BLAS.symv!('U', 1.0, X, x, 0.0, y)
		BLAS.gemv!('T', 1.0, U, y, 0.0, w)
		BLAS.gemv!('N', -1.0, U, w, 1.0, y)
	end
	D = LinearMap{Float64}(custom_mul, n; ismutating=true, issymmetric=true)
	@time (λ_rem, v_rem, nconv, niter, nmult, resid) = eigs(D, nev=1, ncv=20, ritzvec=true, which=:LR, tol=1e-1, v0=v_rem - U*(U'*v_rem))
	@show λ_rem
	@show nmult

	R = [XZ XQ]*V_ - U*Diagonal(λ)
	@show residual = sqrt(2*norm(R)^2 + (n - size(U, 2))*max(λ_rem[1], 0)^2)

	return Xπ, [U Ub], [λ; λb], v_rem[:, 1]
end


QOCS.sdcone_lanczos!(copy(X[:]), cone)
@show cone.subspace_dimension
X = X + Symmetric(0.0001*randn(n, n))
(Xπ, U, λ, v_rem) = my_projection(X, cone.subspace[:, 1:cone.subspace_dimension], cone.λ, randn(size(X, 1)))
X = X + Symmetric(0.0001*randn(n, n))
@time (Xπ, U, λ, v_rem) = my_projection(X, cone.subspace[:, 1:cone.subspace_dimension], cone.λ, v_rem)
X = X + Symmetric(0.1*randn(n, n))
@time (Xπ, U, λ, v_rem) = my_projection(X, cone.subspace[:, 1:cone.subspace_dimension], cone.λ, v_rem)
# @btime my_projection(X, cone.subspace[:, 1:cone.subspace_dimension], cone.λ, v_rem)
# @btime orig_projection(X, cone.subspace[:, 1:cone.subspace_dimension])
cone.subspace_dimension = -2
x2=X[:]
@btime QOCS.sdcone_lanczos!(x2, cone)

X = X + Symmetric(0.0001*randn(n, n))
x1 = X[:]
using Profile
Profile.@profile QOCS.sdcone_lanczos!(x1, cone)
open("prof.txt", "w") do s
	Profile.print(IOContext(s, :displaysize => (24, 500)))
end
#=
cone.subspace_dimension = -1
Profile.@profile QOCS.sdcone_lanczos!(X[:], cone)
open("prof_exact.txt", "w") do s
	Profile.print(IOContext(s, :displaysize => (24, 600)))
end
=#