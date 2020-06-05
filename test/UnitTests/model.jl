
using COSMO, Test, SparseArrays, Random, LinearAlgebra

# This test is precision agnostic and will be run with TestFloat precision
if @isdefined UnitTestFloat
    TestFloat = UnitTestFloat
else
    TestFloat = Float64
end

TestFloat = Float32

rng = Random.MersenneTwister(1872381)

P_int = 4
q_int = 2

P_Float = TestFloat(4.)
q_Float = TestFloat(2.)

P_UInt = UInt64(4)
q_UInt = UInt64(2)

P_mat = rand(rng, TestFloat, 10, 10)
q_mat = rand(rng, TestFloat, 10, 1)


P_sparse = sprand(rng, TestFloat, 10, 10, 0.4)
q_sparse = rand(rng, TestFloat, 10, 1)

A_vec = TestFloat[1 2 3 4]
b_vec = TestFloat(1)


scalar_c = COSMO.Constraint(TestFloat(1.), TestFloat(2.), COSMO.ZeroSet)
constr = COSMO.Constraint(diagm(0 => ones(TestFloat, 10)), rand(rng, TestFloat, 10), COSMO.ZeroSet)


@testset "Model" begin
    # default case
    model = COSMO.Model()
    @test typeof(model) == COSMO.Model{Float64}


    # If integers are used, we convert to model parameter type
    model = COSMO.Model{TestFloat}()
    @test assemble!(model, P_int, q_int, [scalar_c]) == nothing

    model = COSMO.Model{TestFloat}()
    @test assemble!(model, P_Float, q_Float, [scalar_c]) == nothing

    model = COSMO.Model{TestFloat}()
    @test assemble!(model, P_Float, q_Float, [scalar_c]) == nothing


    model = COSMO.Model{TestFloat}()
    @test assemble!(model, P_mat, q_mat, [constr]) == nothing

    model = COSMO.Model{TestFloat}()
    @test assemble!(model, P_sparse, q_sparse, [constr]) == nothing


    # Prevent mixing of precisions
    if TestFloat == Float32
        model = COSMO.Model()
        @test_throws ArgumentError assemble!(model, P_sparse, q_sparse, [constr])

        # test that settings and model have consistent types
        model = COSMO.Model{Float64}()
        settings = COSMO.Settings{Float32}()
        @test_throws ArgumentError assemble!(model, rand(Float64, 10, 10), rand(Float64, 10), [constr], settings = settings)
    end

end
nothing
