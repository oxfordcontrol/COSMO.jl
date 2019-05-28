# Unit Test for by default supported convex sets and their functions
using COSMO, Test, Random, LinearAlgebra

rng = Random.MersenneTwister(13131)

@testset "Convex Sets" begin
tol = 1e-4

    @testset "in_cone Functions" begin
        in_cone_tol = 1e-8
        @test COSMO.in_cone([-1e-6; 0; 1e-6], COSMO.ExponentialCone(), in_cone_tol)
        @test !COSMO.in_cone([-1e-6; 0; -1e-6], COSMO.ExponentialCone(), in_cone_tol)
        x = 2.
        y = -2.
        z = y * exp(x/y)
        @test !COSMO.in_cone([x; y; z], COSMO.ExponentialCone(), in_cone_tol)
        y = 4
        z = y * exp(x/y)
        @test COSMO.in_cone([x; y; z], COSMO.ExponentialCone(), in_cone_tol)

        @test COSMO.in_cone([1e-6; 0; 0], COSMO.PowerCone(0.9), in_cone_tol)
        @test !COSMO.in_cone([-1e-6; 0; 0], COSMO.PowerCone(0.9), in_cone_tol)
        @test !COSMO.in_cone([-1.; 1; 1], COSMO.PowerCone(0.5), in_cone_tol)
        @test COSMO.in_cone([2; 4; 0.1], COSMO.PowerCone(0.5), in_cone_tol)
    end

    @testset "Create and project" begin

    # Zero Cone
    zset = COSMO.ZeroSet(10)
    x = randn(rng,10)
    COSMO.project!(view(x, 1:length(x)), zset)
    @test norm(x,Inf) == 0.

    # Positive Orthant R+
    nonnegatives = COSMO.Nonnegatives(10)
    x = randn(rng, 10)
    COSMO.project!(view(x, 1:length(x)), nonnegatives)
    @test minimum(x) >= 0.

    # Box
    l = -1 * ones(10)
    u = 1 * ones(10)
    box = COSMO.Box(l, u)
    box = COSMO.Box{Float64}(10)
    box.l .= l
    box.u .= u
    x = 100 * randn(rng, 10)
    COSMO.project!(view(x, 1:length(x)), box)
    @test minimum(x) >= -1. && maximum(x) <= 1.

    # Second Order (Lorentz) cones
    soc = COSMO.SecondOrderCone(10)
    x = 10 * randn(rng, 9)
    t = norm(x, 2) - 0.5
    x = [t; x]

    COSMO.project!(view(x, 1:length(x)), soc)
    @test norm(x[2:10], 2) <= x[1]

    # Positive Semidefinite cones
    psd = COSMO.PsdCone(16)
    X = Symmetric(randn(rng, 4, 4))
    X = X - 4 * Matrix(1.0I, 4, 4)
    x = vec(X)
    COSMO.project!(view(x, 1:length(x)), psd)
    @test minimum(eigen(reshape(x, 4, 4)).values) >= -1e-9

    C = COSMO.CompositeConvexSet([COSMO.ZeroSet(10), COSMO.Nonnegatives(10)])
    x = -rand(20)
    xs = COSMO.SplitVector(x, C)
    COSMO.project!(xs, C)
    @test norm(x, Inf) == 0.

    # Exponential cones
    exp_cone = COSMO.ExponentialCone()
    exp_test_passed = true
    for i = 1:100
        x = -25 .+ rand(rng, 3) .* 50
        x_in = copy(x)
        COSMO.project!(x, exp_cone)
        !COSMO.in_cone(x, exp_cone, 1e-4) && (exp_test_passed = false)
    end
    @test exp_test_passed

    # 3-d Power cones
    pow_test_passed = true
    for i = 1:100
        pow_cone = COSMO.PowerCone(0.1 .+ rand(rng) * 0.85)
        x = -25 .+ rand(rng, 3) .* 50
        x_in = copy(x)
        COSMO.project!(x, pow_cone)
        if !COSMO.in_cone(x, pow_cone, 1e-4)
            @show(pow_cone.α)
            @show(x, x_in)
        end

        !COSMO.in_cone(x, pow_cone, 1e-4) && (pow_test_passed = false)
    end
    @test pow_test_passed
    end


    @testset "in_dual Functions" begin

    # Dual of zero cone
    x = randn(rng, 10)
    convex_set = COSMO.ZeroSet(10)
    @test COSMO.in_dual(view(x, 1:length(x)), convex_set, tol)

    # Dual of Positive Orthant R+ (self-dual)
    xpos = rand(rng, 10)
    xneg = -rand(rng, 10)
    xzeros = zeros(10)
    convex_set = COSMO.Nonnegatives(10)
    @test COSMO.in_dual(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_dual(view(xneg,1:length(xneg)), convex_set, tol)
    @test COSMO.in_dual(view(xzeros,1:length(xzeros)), convex_set,tol)

    #TODO: Dual of Box [important!]

    # Dual of Second Order Cone (self-dual)
    tol = 1e-4
    x = randn(rng, 9)
    t = norm(x, 2)
    xpos = [t + 0.5; x]
    xneg = [t - 0.5; x]
    convex_set = COSMO.SecondOrderCone(10)
    @test COSMO.in_dual(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_dual(view(xneg, 1:length(xneg)), convex_set, tol)

    # Dual of Positive Semidefinite Cone (Square) (self-dual)
    tol = 1e-4
    X = Symmetric(randn(rng, 4, 4))
    Xpos = X + 5 * Matrix(1.0I, 4, 4)
    Xneg = X - 5 * Matrix(1.0I, 4, 4)
    xpos = vec(Xpos)
    xneg = vec(Xneg)
    convex_set = COSMO.PsdCone(16)
    @test COSMO.in_dual(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_dual(view(xneg, 1:length(xneg)), convex_set, tol)

    # Dual of Positive Semidefinite Cone (Triangle)
    convex_set = COSMO.PsdConeTriangle(10)
    xpos = zeros(10)
    xneg = zeros(10)
    COSMO.extract_upper_triangle!(Xpos, xpos, sqrt(2))
    COSMO.extract_upper_triangle!(Xneg, xneg, sqrt(2))
    @test COSMO.in_dual(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_dual(view(xneg, 1:length(xneg)), convex_set, tol)

    # Dual of Exponential Cone
    convex_set = COSMO.ExponentialCone()
    xpos = [0; rand(rng, 2)]
    x = -2.
    y = 3.
    z = -x * exp(y / x) / exp(1) + 0.1
    xpos2 = [x; y; z]
    xneg = rand(rng, 3)
    xneg2 = [x; y; z - 0.2]
    @test COSMO.in_dual(xpos, convex_set, tol)
    @test COSMO.in_dual(xpos2, convex_set, tol)
    @test !COSMO.in_dual(xneg, convex_set, tol)
    @test !COSMO.in_dual(xneg2, convex_set, tol)

    # Dual of Power Cone
    convex_set = COSMO.PowerCone(0.5)
    xpos = [1.; 2.; 2.5]
    xneg = [1.; 2.; 3]
    @test COSMO.in_dual(xpos, convex_set, tol)
    @test !COSMO.in_dual(xneg, convex_set, tol)
    end

    @testset "in_pol_recc Functions" begin

    # Polar Recession cone of zero cone
    xpos = zeros(10)
    xneg = randn(rng, 10)
    convex_set = COSMO.ZeroSet(10)
    @test COSMO.in_pol_recc(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_pol_recc(view(xneg, 1:length(xneg)), convex_set, tol)

    # Polar Recession cone of Positive Orthant R+
    xpos = -rand(rng, 10)
    xneg = rand(rng, 10)
    convex_set = COSMO.Nonnegatives(10)
    @test COSMO.in_pol_recc(view(xpos,1:length(xpos)), convex_set, tol)
    @test !COSMO.in_pol_recc(view(xneg,1:length(xneg)), convex_set, tol)

    # Polar Recc of Second Order Cone
    tol = 1e-4
    x = randn(rng, 9)
    t = norm(x, 2)
    xpos = [-t - 0.5; x]
    xneg = [-t + 0.5; x]
    convex_set = COSMO.SecondOrderCone(10)
    @test COSMO.in_pol_recc(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_pol_recc(view(xneg, 1:length(xneg)), convex_set, tol)

    # Polar Recc of Positive Semidefinite Cone (Square)
    tol = 1e-4
    X = Symmetric(randn(rng, 4, 4))
    Xpos = X - 20 * Matrix(1.0I, 4, 4)
    Xneg = X + 4 * Matrix(1.0I, 4, 4)
    xpos = vec(Xpos)
    xneg = vec(Xneg)
    convex_set = COSMO.PsdCone(16)
    @test COSMO.in_pol_recc(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_pol_recc(view(xneg, 1:length(xneg)), convex_set, tol)

    # Polar Recc of Positive Semidefinite Cone (Triangle)
    COSMO.extract_upper_triangle!(Xpos, xpos, sqrt(2))
    COSMO.extract_upper_triangle!(Xneg, xneg, sqrt(2))
    convex_set = COSMO.PsdConeTriangle(10)
    @test COSMO.in_pol_recc(view(xpos, 1:length(xpos)), convex_set, tol)
    @test !COSMO.in_pol_recc(view(xneg, 1:length(xneg)), convex_set, tol)

    # Polar Recc of Exponential Cone
    convex_set = COSMO.ExponentialCone()
    xpos = [0; -rand(rng, 2)]
    x = 2.
    y = -3.
    z = -0.2641699972477976
    xpos2 = [x; y; z]
    xneg = -rand(rng, 3)
    xneg2 = [x; y; -z]
    @test COSMO.in_pol_recc(xpos, convex_set, tol)
    @test COSMO.in_pol_recc(xpos2, convex_set, tol)
    @test !COSMO.in_pol_recc(xneg, convex_set, tol)
    @test !COSMO.in_pol_recc(xneg2, convex_set, tol)

    # Polar Recc of power Cone
    @test COSMO.in_pol_recc([-1; -2; -2.5], COSMO.PowerCone(0.5), tol)
   end

   @testset "Set Utilities" begin
    convex_set = COSMO.PsdConeTriangle(3)
    composite_set = COSMO.CompositeConvexSet([COSMO.ZeroSet(2), COSMO.Nonnegatives(1)])
    @test COSMO.get_subset(convex_set, 1) == convex_set
    @test_throws DimensionMismatch COSMO.get_subset(convex_set, 2)
    @test COSMO.get_subset(composite_set, 1) == COSMO.ZeroSet(2)
    @test_throws ErrorException COSMO.CompositeConvexSet([convex_set; composite_set])
   end

end


nothing
