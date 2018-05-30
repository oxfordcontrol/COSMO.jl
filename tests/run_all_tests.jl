using Base.Test

@testset "All tests" begin
    include("./UnitTests/InfeasibilityTests/runTests.jl")
    include("./UnitTests/ScalingTests/runTests.jl")
end;