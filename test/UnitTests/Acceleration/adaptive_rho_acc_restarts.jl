# Test whether the accelerator is restarted when rho is adapted (the ADMM operator changed)

using COSMO, SparseArrays, LinearAlgebra, Test


@testset "Rho adaption triggers restarts" begin
    q = [1; 1.];
    P = sparse([4. 1; 1 2]);
    A = [1. 1; 1 0; 0 1];
    l = [1.; 0; 0];
    u = [1; 0.7; 0.7];

    Aa = [-A; A]
    ba = [u; -l]
    constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);

    accelerator = with_options(AndersonAccelerator{Float64, NoRegularizer, Type2, RollingMemory}, mem = 5, safeguarded = false, activation_reason = ImmediateActivation())
    settings = COSMO.Settings(adaptive_rho_interval = 23, rho = 1e-4, eps_abs = 1e-5, eps_rel = 1e-5, accelerator = accelerator);
    model = COSMO.Model();
    assemble!(model, P, q, constraint1, settings = settings);
    res = COSMO.optimize!(model);
    
    @test length(filter(x -> x[2] == :rho_adapted, model.accelerator.restart_iter)) == num_rho_adaptions(model)
end