using LinearAlgebra
using Random
using Test
using BenchmarkTools
using Profile

function generate_test_matrix(n, n_pos, n_small, gap)
    # Generates the n x n symmetric matrix A with
    # n_pos eigenvalues logarithmically on [0, 1]
    # n_small eigenvalues uniformly on [-gap, gap]
    # and the rest uniformly on [-1, 0].
    n_neg = n - n_pos - n_small;
    @assert n_neg >=0 "Size of matrix is too small"

    eigenvalues_positive = 10.0.^range(-10, 0, length=n_pos)
    eigenvalues_negative = range(-1, 0, length=n_neg);
    eigenvalues_small = gap*range(-1, 1, length=n_small);

    Q = Matrix(qr(randn(n, n)).Q)
    l = [eigenvalues_positive; eigenvalues_negative; eigenvalues_small]
    return Symmetric(Q'*Diagonal(l)*Q)
end

Random.seed!(1)
n = 2000; n_pos = 3; n_small = 2;
A = generate_test_matrix(n, n_pos, n_small, 1e-10);
X0 = randn(n, n_pos)
max_iter = 1000
tol=1e-7
verbosity=1
function check_residuals(A, data, barrier = 0.0)
    data.X = data.X[:, data.λ .>= barrier]
    data.λ = data.λ[data.λ .>= barrier]
    # The code does not really guarrantie the requested tolerance. See how active_indices are defined in the Code.
    @test all(COSMO.column_norms(A*data.X - data.X*Diagonal(data.λ)) .< 4*tol)
    @test size(data.X, 2) >= sum(eigvals(A) .> tol)
end

@testset "LOBPCG - Basic Tests" begin
    @testset "Matrix with predefined eigenvalues" begin
        data, flag = COSMO.lobpcg(A, copy(X0), max_iter, tol=tol,
            λ_barrier=0.0, verbosity=verbosity)
        check_residuals(A, data, 0.0)
        # Test which=:smallest
        data, flag = COSMO.lobpcg(-A, copy(X0), max_iter, tol=tol,
            λ_barrier=0.0, verbosity=verbosity, which=:smallest)
        data.λ = -data.λ
        check_residuals(A, data, 0.0)
    end
    @testset "Excessive desired eigenvalues" begin
        data, flag = COSMO.lobpcg(-A, copy(X0), max_iter, tol=tol,
            λ_barrier=0.0, verbosity=verbosity)
        @test flag == :excessive_column_num
    end

    @testset "Matrix with no positive eigenvalues" begin
        A = randn(n, 50); A = Symmetric(-A*A');
        data, flag = COSMO.lobpcg(A, copy(X0), max_iter, tol=tol,
            λ_barrier=0.0,
            verbosity=verbosity)
        check_residuals(A, data, 0.0)
    end

    @testset "Matrix of all zeros" begin
        A = zeros(n, n)
        data, flag = COSMO.lobpcg(A, copy(X0), max_iter, tol=tol,
            λ_barrier=0.0, verbosity=verbosity)
        check_residuals(A, data, 0.0)
    end
end