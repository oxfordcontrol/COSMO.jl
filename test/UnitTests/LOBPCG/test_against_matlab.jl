using MATLAB
using LinearAlgebra, Random
using Test
include("../../../src/lobpcg.jl")

Random.seed!(1)
n = 2500 # Matrix size
k = 100  # Number of deisred eigenvalues
tol = 1e-4
A = Matrix(Symmetric(randn(n, n)))
X0 = randn(n, k)
max_iter = 10
verbosity = 2

function compute_eigenvector_error(X1, X2)
    error_sum = 0
    for i = 1:k
        error_sum += min(norm(X1[:, i] - X2[:, i]), norm(X1[:, i] + X2[:, i]))
    end
    return error_sum
end

@testset "Random Matrix: Compare with MATLAB's implementation" begin
    data, converged = lobpcg(A, copy(X0), tol=tol, max_iter, which=:smallest, verbosity=verbosity)
    blockVectorX, lambda, failureFlag = mxcall(:lobpcg_stripped, 3, X0, A, tol, max_iter)
    # show(stdout, "text/plain", data.X); println()
    # show(stdout, "text/plain", blockVectorX); println()
    @test norm(lambda - data.λ[end:-1:1]) < 1e-8
    @test compute_eigenvector_error(data.X[:, end:-1:1], blockVectorX) < 1e-8
end

nothing