using COSMO, Random, Test, Pkg
rng = Random.MersenneTwister(12345)

@testset "Updatable Q" begin
    n = 700
    k = 200

    X = Symmetric(randn(rng, n, n));
    A = randn(n, k)
    F = UpdatableQ(A)
    add_columns!(F, X*F.Q1)
    @test F.Q1'*F.Q1 ≈ I
    R = F.Q1'*[A X*A]
    R1 = R[1:size(R,2), :]
    R2 = R[size(R,2)+1:end, :]
    @test UpperTriangular(R1) ≈ R1
    show(stdout, "text/plain", R2); println()
    @test norm(R2) < 1e-9
end
