using Random, Test, LinearAlgebra

verbosity = 0 # Change to 2 to see LOBPCG's output
Random.seed!(1)
λ = [1:6; -50.00:-0.00]
n = length(λ)
Q = Matrix(qr(randn(n, n)).Q)
X = Q*Diagonal(λ)*Q'

n_ = Int(n*(n + 1)/2)

function exact_projection(A::Matrix)
    l, V = eigen(Symmetric(A))
    return Symmetric(V*Diagonal(max.(l, 0))*V')
end

function lobpcg_unit_test!(cone, X::Matrix)
    cone.exact_projections = 0
    cone.lobpcg_iterations = 0
    x = COSMO.extract_upper_triangle(X)
    X_proj_exact = exact_projection(X)
    COSMO.project!(x, cone)
    X_proj = COSMO.populate_upper_triangle(x)

    @test norm(X_proj - X_proj_exact) < 1e-5 # This is tolerance is hand-wavy
end

function set_tolerance!(cone, tol)
    # Currently, the tolerance of LOBPCG is determined as a
    # monotonically decreasing function of LOBPCG's iteration number.
    # This function (i.e. set_tolerance!) uses binary search to find
    # the iteration number that gives a close tolerance to the specified tol.
    iter_min = 0
    iter_max = 1e12
    for i = 1:100
        cone.iteration = Int(round((iter_min + iter_max)/2))
        if COSMO.get_tolerance(cone) < tol
            iter_max = cone.iteration
        else
            iter_min = cone.iteration
        end
    end
    # @show COSMO.get_tolerance(cone), tol
end

function lobpcg_unit_tests(X::Matrix)
    cone = COSMO.PsdConeTriangleLOBPCG(n_, verbosity=verbosity, max_iter=200)
    @testset "Initial call (exact eigendecomposition)" begin
        lobpcg_unit_test!(cone, X)
        @test cone.exact_projections == 1
        @test cone.lobpcg_iterations == 0
    end

    @testset "Warm starting - already solved to tolerance" begin
        set_tolerance!(cone, 1e-1)
        lobpcg_unit_test!(cone, X + 1e-6*Symmetric(randn(n, n)))
        @test cone.exact_projections == 0
        @test cone.lobpcg_iterations == 1
    end
    @testset "Warm starting - a few iterations required" begin
        set_tolerance!(cone, 1e-7)
        # Perturb sligtly the matrix under projection
        lobpcg_unit_test!(cone, X + 1e-6*Symmetric(randn(n, n)))
        @test cone.exact_projections == 0
        @test cone.lobpcg_iterations >= 2
        @test cone.lobpcg_iterations < 10
    end

    @testset "Warm starting - expand subspace" begin
        set_tolerance!(cone, 1e-7)
        # Shift the spectrum
        if cone.is_subspace_positive
            lobpcg_unit_test!(cone, X + 3.1*I)
        else
            lobpcg_unit_test!(cone, X - 3.1*I)
        end
        @test size(cone.U, 2) > 6
        @test cone.exact_projections == 0
        @test cone.lobpcg_iterations < 35
    end

    @testset "Warm starting - contract subspace" begin
        set_tolerance!(cone, 1e-7)
        # Shift the spectrum
        if cone.is_subspace_positive
            lobpcg_unit_test!(cone, X - 1.1*I)
        else
            lobpcg_unit_test!(cone, X + 1.1*I)
        end
        @test size(cone.U, 2) < 6
        @test cone.exact_projections == 0
        @test cone.lobpcg_iterations < 35
    end
end

@testset "LOBPCG PSD Projection" begin
    @testset "Track Positive Eigenspace" begin
        lobpcg_unit_tests(X)
    end
    @testset "Track Negative Eigenspace" begin
        lobpcg_unit_tests(-X)
    end
end

nothing
