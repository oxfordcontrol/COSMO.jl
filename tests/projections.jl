# Test file to check the projections
workspace()
include("../Projections.jl")
include("../Helper.jl")
using Base.Test, Helper

rng = MersenneTwister(1234)

nn = 1000


@testset "Test projection onto zero cone" begin
        x = rand(rng,rand(rng,1:100),1)
        @test Projections.zeroCone(x) == zeros(size(x,1),1)
end


@testset "Test projection onto nonNegativeOrthant" begin
    for iii = 1:nn
        dim = rand(rng,1:100)
        x = rand(rng,-1:0.001:1,dim,1)
        @test minimum(Projections.nonNegativeOrthant(x)) >= 0
    end
end

@testset "Test projection onto Box" begin
    for iii = 1:nn
        dim = rand(rng,1:100)
        vMin = -1000
        vMax = 1000
        x = vMin+rand(rng,dim,1)*(vMax-vMin)
        u = (vMin/5)+rand(rng,1e-6:1e-3:1e6,dim,1)*((vMax-vMin)/5)
        l = -u
        xProj = Projections.box(x,l,u)

        projectionFailed = false
        for kkk = 1:dim
            if xProj[kkk] < l[kkk] || xProj[kkk] > u[kkk]
                projectionFailed = true
            end
        end
        @test projectionFailed == false
    end
end

@testset "Test projection onto second-order cone" begin

    for iii = 1:nn
        # create random vector

        dim = rand(rng,1:100)
        vMin = -2
        vMax = 2
        x = vMin+rand(rng,dim,1)*(vMax-vMin)
        t = rand(rng)
        # project matrix onto positive semidefinite cone
        xNew,tNew = Projections.secondOrderCone(x,t)

        # check if relation holds {(t,x) | ||x||_2 <= t}
        @test (norm(xNew,2) - tNew) < 1e-14

    end
end

@testset "Test projection onto pos-def cone" begin

  for iii = 1:nn

    # create random matrix
    dim = rand(rng,1:100)
    A = rand(rng,dim,dim)
    A = full(Symmetric(A))
    # project matrix onto positive semidefinite cone
    a = vec(A)
    X = Projections.sdcone(a,dim)

    # check positive semidefiniteness
    @test isNumericallyPosSemDef(reshape(X,dim,dim),-1e-10) == true

  end

end