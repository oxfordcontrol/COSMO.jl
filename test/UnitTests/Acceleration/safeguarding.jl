# test safeguarding related functions

using COSMO, SparseArrays, LinearAlgebra, Test

@testset "Acceleration safeguarding" begin
    q = [1; 1.];
    P = sparse([4. 1; 1 2]);
    A = [1. 1; 1 0; 0 1];
    l = [1.; 0; 0];
    u = [1; 0.7; 0.7];

    Aa = [-A; A]
    ba = [u; -l]
    constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);
    m, n = size(Aa)
    # Next, we define the settings object, the model and then assemble everything:
    # model = COSMO.Model();
    # assemble!(model, P, q, constraint1, settings = COSMO.Settings(max_iter = 10, scaling = 0));
	# model.rws = COSMO.ResidualWorkspace{Float64}(m, n, model.p.C) 
    # rws = model.rws
    # # make sure rws.ν works as a view
    # @. rws.sol = rand()
    # @test rws.sol[n+1:end] == rws.ν
    
    # # do ten iterations
    # result = COSMO.optimize!(model);
    # w0 = deepcopy(model.vars.w)
    # # do one more iteratios
    # model.settings.max_iter = 1
    # result = COSMO.optimize!(model);
    # w1 = deepcopy(model.vars.w)

    # # now put w0 through fixed-point residual function 
    # fp_res = COSMO.fixed_point_residual_norm(model.rws, model, w0)
    # @test fp_res == norm(w0 - w1, 2)


    # now check algorithm behaviour with safeguarding
    model = COSMO.Model();
    accelerator = with_options(AndersonAccelerator{Float64, NoRegularizer, Type2, RollingMemory}, mem = 5, safeguarded = true, activation_reason = ImmediateActivation())
    settings = COSMO.Settings(eps_abs = 1e-6, eps_rel = 1e-6, accelerator = accelerator)
    assemble!(model, P, q, constraint1, settings = settings);
    result = COSMO.optimize!(model);

end