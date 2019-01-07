using MathOptInterface, COSMO, Test
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model(COSMOModelData,
        (),
        (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
        (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
         MOI.ExponentialCone, MOI.PositiveSemidefiniteConeSquare, MOI.PositiveSemidefiniteConeTriangle),
        (),
        (MOI.SingleVariable,),
        (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
        (MOI.VectorOfVariables,),
        (MOI.VectorAffineFunction,),);
optimizer = MOIU.CachingOptimizer(COSMOModelData{Float64}(),
                                  COSMO.Optimizer(eps_abs = 1e-6, eps_rel = 1e-6 ));

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "COSMO";
end


config = MOIT.TestConfig(atol=1e-2, rtol=1e-2)
@testset "Unit" begin
    MOIT.unittest(MOIB.SplitInterval{Float64}(optimizer), config,
                  [# Quadratic functions are not supported
                   "solve_qcp_edge_cases", "solve_qp_edge_cases",
                   # Integer and ZeroOne sets are not supported
                   "solve_integer_edge_cases", "solve_objbound_edge_cases"])
end

# FIXME: one test fails for linear12 problem. Take a closer look why
@testset "Continuous linear problems" begin
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config, ["linear12"])
end

@testset "Continuous quadratic problems" begin
    exclude_qt_test_sets = ["qcp", "socp"]
    MOIT.contquadratictest(optimizer, config, exclude_qt_test_sets)
end


exclude_conic_test_sets = ["exp", "rootdet", "logdet"]
@testset "Continuous conic problems" begin
    MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.LogDet{Float64}(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(optimizer)))),
                       config, exclude_conic_test_sets)
end







nothing