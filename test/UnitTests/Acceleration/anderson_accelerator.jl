# Test the different anderson accelerator types

using COSMO, SparseArrays, LinearAlgebra, Test

anderson_types = [
    AndersonAccelerator{Float64, Type1, RollingMemory, NoRegularizer}; 
    AndersonAccelerator{Float64, Type1, RestartedMemory, NoRegularizer}; 
    AndersonAccelerator{Float64, Type1, RollingMemory, TikonovRegularizer}; 
    AndersonAccelerator{Float64, Type1, RestartedMemory, TikonovRegularizer}; 
    AndersonAccelerator{Float64, Type1, RollingMemory, FrobeniusNormRegularizer}; 
    AndersonAccelerator{Float64, Type1, RestartedMemory, FrobeniusNormRegularizer}; 
    AndersonAccelerator{Float64, Type2{NormalEquations}, RollingMemory, NoRegularizer}; 
    AndersonAccelerator{Float64, Type2{NormalEquations}, RestartedMemory, NoRegularizer}; 
    AndersonAccelerator{Float64, Type2{NormalEquations}, RollingMemory, TikonovRegularizer}; 
    AndersonAccelerator{Float64, Type2{NormalEquations}, RestartedMemory, TikonovRegularizer}; 
    AndersonAccelerator{Float64, Type2{NormalEquations}, RollingMemory, FrobeniusNormRegularizer}; 
    AndersonAccelerator{Float64, Type2{NormalEquations}, RestartedMemory, FrobeniusNormRegularizer}; 
    AndersonAccelerator{Float64, Type2{QRDecomp}, RestartedMemory, NoRegularizer}
] 

@testset "AndersonAccelerator" begin

    @testset "Types" begin
        q = [1; 1.];
        P = sparse([4. 1; 1 2]);
        A = [1. 1; 1 0; 0 1];
        l = [1.; 0; 0];
        u = [1; 0.7; 0.7];

        Aa = [-A; A]
        ba = [u; -l]
        constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);

        for at in anderson_types
            settings = COSMO.Settings(accelerator = at);
            model = COSMO.Model();
            assemble!(model, P, q, constraint1, settings = settings);
            res = COSMO.optimize!(model);
            @test res.status == :Solved
        end
    end
end