# Unit Test for by default supported convex sets and their functions
using COSMO, Test, Random,LinearAlgebra

rng = Random.MersenneTwister(13131)

@testset "Convex Sets" begin
tol = 1e-4


    @testset "Create and project" begin

    # Zero Cone
    zset = COSMO.ZeroSet(10)
    x = randn(rng,10)
    COSMO.project!(view(x,:),zset)
    @test norm(x,Inf) == 0.

    # Positive Orthant R+
    nonnegatives = COSMO.Nonnegatives(10)
    x = randn(rng,10)
    COSMO.project!(view(x,:),nonnegatives)
    @test minimum(x) >= 0.

    # Box
    l = -1*ones(10)
    u = 1*ones(10)
    box = COSMO.Box(l,u)
    x = 100*randn(rng,10)
    COSMO.project!(view(x,:),box)
    @test minimum(x) >= -1. && maximum(x) <= 1.

    # Second Order (Lorentz) cones
    soc = COSMO.SecondOrderCone(10)
    x = 10*randn(rng,9)
    t = norm(x,2) - 0.5
    x = [t;x]

    COSMO.project!(view(x,:),soc)
    @test norm(x[2:10],2) <= x[1]

    # Positive Semidefinite cones
    psd = COSMO.PsdCone(16)
    X = randn(rng,4,4)
    X = X*X' - 4*Matrix(1.0I,4,4)
    x = vec(X)
    COSMO.project!(view(x,:),psd)
    @test minimum(eigen(reshape(x,4,4)).values) >= -1e-9
    end


    @testset "indual Functions" begin

    # Dual of zero cone
    x = randn(rng,10)
    convexSet = COSMO.ZeroSet(10)
    @test COSMO.indual(view(x,:),convexSet,tol)

    # Dual of Positive Orthant R+ (self-dual)
    xpos = rand(rng,10)
    xneg = -rand(rng,10)
    xzeros = zeros(10)
    convexSet = COSMO.Nonnegatives(10)
    @test COSMO.indual(view(xpos,:),convexSet,tol)
    @test !COSMO.indual(view(xneg,:),convexSet,tol)
    @test COSMO.indual(view(xzeros,:),convexSet,tol)

    #TODO: Dual of Box [important!]

    # Dual of Second Order Cone (self-dual)
    tol = 1e-4
    x = randn(rng,9)
    t = norm(x,2)
    xpos = [t+0.5;x]
    xneg = [t-0.5;x]
    convexSet = COSMO.SecondOrderCone(10)
    @test COSMO.indual(view(xpos,:),convexSet,tol)
    @test !COSMO.indual(view(xneg,:),convexSet,tol)



    # Dual of Positive Semidefinite Cone (self-dual)
    tol = 1e-4
    X = randn(rng,4,4)
    Xpos = X*X' + 4*Matrix(1.0I,4,4)
    Xneg = X*X' - 4*Matrix(1.0I,4,4)
    xpos = vec(Xpos)
    xneg = vec(Xneg)
    convexSet = COSMO.PsdCone(16)
    @test COSMO.indual(view(xpos,:),convexSet,tol)
    @test !COSMO.indual(view(xneg,:),convexSet,tol)

    end

    @testset "inPolRec Functions" begin

    # Polar Recession cone of zero cone
    xpos = zeros(10)
    xneg = randn(rng,10)
    convexSet = COSMO.ZeroSet(10)
    @test COSMO.inrecc(view(xpos,:),convexSet,tol)
    @test !COSMO.inrecc(view(xneg,:),convexSet,tol)

    # Polar Recession cone of Positive Orthant R+
    xpos = -rand(rng,10)
    xneg = rand(rng,10)
    convexSet = COSMO.Nonnegatives(10)
    @test COSMO.inrecc(view(xpos,:),convexSet,tol)
    @test !COSMO.inrecc(view(xneg,:),convexSet,tol)

    #TODO: Polar Recc of Box [important!]

    # Polar Recc of Second Order Cone
    tol = 1e-4
    x = randn(rng,9)
    t = norm(x,2)
    xpos = [-t-0.5;x]
    xneg = [-t+0.5;x]
    convexSet = COSMO.SecondOrderCone(10)
    @test COSMO.inrecc(view(xpos,:),convexSet,tol)
    @test !COSMO.inrecc(view(xneg,:),convexSet,tol)



    # Polar Recc of Positive Semidefinite Cone
    tol = 1e-4
    X = randn(rng,4,4)
    Xpos = X*X' - 20*Matrix(1.0I,4,4)
    Xneg = X*X' + 4*Matrix(1.0I,4,4)
    xpos = vec(Xpos)
    xneg = vec(Xneg)
    convexSet = COSMO.PsdCone(16)
    @test COSMO.inrecc(view(xpos,:),convexSet,tol)
    @test !COSMO.inrecc(view(xneg,:),convexSet,tol)

    end

end


nothing
