# Test file to check the projections
# workspace()
# include("../src/Solver.jl")
using Test, COSMOTestUtils, COSMO.Projections, LinearAlgebra, SparseArrays

rng = Random.MersenneTwister(1234)

nn = 1


@testset "Test projection onto zero cone" begin
        x = rand(rng, rand(rng, 1:100))
        Projections.zeroCone!(x, 1, length(x))
        @test x == zeros(size(x, 1))
end


@testset "Test projection onto nonNegativeOrthant" begin
    for iii = 1:nn
        dim = rand(rng, 1:100)
        x = rand(rng, -1:0.001:1, dim)
        Projections.nonNegativeOrthant!(x, 1, length(x))
        @test minimum(x) >= 0
    end
end

@testset "Test projection onto Box" begin
    for iii = 1:nn
        dim = rand(rng, 1:100)
        vMin = -1000
        vMax = 1000
        x = vMin + rand(rng, dim, 1) * (vMax - vMin)
        u = (vMin / 5) + rand(rng, 1e-6:1e-3:1e6, dim ,1) * ((vMax - vMin) / 5)
        l = -u
        xProj = Projections.box(x, l, u)

        projection_failed = false
        for kkk = 1:dim
            if xProj[kkk] < l[kkk] || xProj[kkk] > u[kkk]
                projection_failed = true
            end
        end
        @test projection_failed == false
    end
end

@testset "Test projection onto second-order cone" begin

    for iii = 1:nn
        # create random vector

        dim = rand(rng, 1:100)
        vMin = -2
        vMax = 2
        x = vMin + rand(rng, dim, 1) * (vMax - vMin)
        t = rand(rng)
        x = [x; t]
        # project matrix onto positive semidefinite cone
        xNew = Projections.secondOrderCone(x)

        # check if relation holds {(t,x) | ||x||_2 <= t}
        @test (norm(xNew[2:end], 2) - xNew[1]) < 1e-14

    end
end

@testset "Test projection onto pos-def cone" begin

  for iii = 1:nn

    # create random matrix
    dim = rand(rng, 1:100)
    A = rand(rng, dim, dim)
    A = full(Symmetric(A))
    # project matrix onto positive semidefinite cone
    a = vec(A)
    X = Projections.sdcone(a)

    # check positive semidefiniteness
    @test minimum(eig(reshape(X, dim, dim))[1]) >= -1e-10

  end

end
nothing
