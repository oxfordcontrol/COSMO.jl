# Test that the settings adaptive_rho_max_adaptions is correctly implemented


using COSMO, SparseArrays, LinearAlgebra, Test


@testset "Max number of rho adaptions" begin
    q = [1; 1.];
    P = sparse([4. 1; 1 2]);
    A = [1. 1; 1 0; 0 1];
    l = [1.; 0; 0];
    u = [1; 0.7; 0.7];

    Aa = [-A; A]
    ba = [u; -l]
    constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);

    # Next, we define the settings object, the model and then assemble everything:
    settings = COSMO.Settings(verbose = false, adaptive_rho_interval = 25, adaptive_rho_max_adaptions = 2, rho = 1e-6, eps_abs = 1e-6, eps_rel = 1e-4);
    model = COSMO.Model();
    assemble!(model, P, q, constraint1, settings = settings);
    res = COSMO.optimize!(model);

    @test num_rho_adaptions(model) == 2

    # Solve again but restrict number of rho updates to 1
    settings = COSMO.Settings(verbose = false, adaptive_rho_interval = 25, adaptive_rho_max_adaptions = 1, rho = 1e-6, eps_abs = 1e-4, eps_rel = 1e-4);
    model = COSMO.Model();
    assemble!(model, P, q, constraint1, settings = settings);
    res = COSMO.optimize!(model);

    @test num_rho_adaptions(model) == 1
end