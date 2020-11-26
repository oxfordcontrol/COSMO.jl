# Test whether the the delayed start of the accelerator works
using COSMO, SparseArrays, LinearAlgebra, Test


@testset "Delayed accelerator start" begin
    
    q = [1; 1.];
    P = sparse([4. 1; 1 2]);
    A = [1. 1; 1 0; 0 1];
    l = [1.; 0; 0];
    u = [1; 0.7; 0.7];

    Aa = [-A; A]
    ba = [u; -l]
    constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);
    
    # acceleration delayed by certain number of iterations
    delayed_iter = 30
    accelerator = with_options(AndersonAccelerator{Float64, NoRegularizer, Type2, RollingMemory}, mem = 5, safeguarded = false, activation_reason = IterActivation(delayed_iter))
    settings = COSMO.Settings(eps_abs = 1e-5, eps_rel = 1e-5, accelerator = accelerator);
    model = COSMO.Model();
    assemble!(model, P, q, constraint1, settings = settings);
    res = COSMO.optimize!(model);
    
    @test delayed_iter == model.accelerator.acceleration_status[1][1]
   
    
    # acceleration delayed by accuracy of the algorithm 
    check_termination = 40
    delayed_acc = 1e-1
    accelerator = with_options(AndersonAccelerator{Float64, NoRegularizer, Type2, RollingMemory}, mem = 5, safeguarded = false, activation_reason = AccuracyActivation(delayed_acc))
    settings = COSMO.Settings(check_termination = check_termination, eps_abs = 1e-5, eps_rel = 1e-5, accelerator = accelerator);
    model = COSMO.Model();
    assemble!(model, P, q, constraint1, settings = settings);
    res = COSMO.optimize!(model);
    
    @test model.accelerator.acceleration_status[1][1] == check_termination + 1 # first time the accuracy is checked
    
end