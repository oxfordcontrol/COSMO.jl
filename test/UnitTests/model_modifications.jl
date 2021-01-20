# This tests the behaviour of multiple successive solves with problem data updates

using Test, LinearAlgebra, SparseArrays
using COSMO

# This test is precision agnostic and will be run with T precision
if !@isdefined(UnitTestFloats)
    UnitTestFloats = [Float64] #if not run in full test setup, just do it for one float type
end

T = Float64

@testset "Model modifications" begin

# for T in UnitTestFloats
    # calling optimize!() twice should yield the same results (with less iterations)
    q = [1; 1.];
    P = sparse([4. 1; 1 2]);
    A = [1. 1; 1 0; 0 1];
    l = [1.; 0; 0];
    u = [1; 0.7; 0.7];
    Aa = [-A; A]
    ba = [u; -l]
    constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);
    model = COSMO.Model();
    assemble!(model, P, q, constraint1, settings= COSMO.Settings(check_termination = 1));
    @test model.states.IS_ASSEMBLED = true
    res = COSMO.optimize!(model);
    res2 = COSMO.optimize!(model);
    @test abs(res.obj_val - res2.obj_val) <= 1e-3
    @test res2.iter <= res.iter

    # updating q should solve correctly without reassembling the KKT matrix
    model = COSMO.Model();
    assemble!(model, P, q, constraint1);
    res = COSMO.optimize!(model);
    # update q and compare to known solution
    q_new = [2.; 3]
    COSMO.update!(model, q = q_new)
    res = COSMO.optimize!(model);
    @test abs(res.obj_val - 3.5) < 1e-3
    @test norm(res.x - [0.5; 0.5], 2) < 1e-3

    # updating b should solve correctly without having to reassemble the KKT matrix, i.e. same factorisation time
    # min   x_1 + x_2
    # s.t.  x_1 >= 2 --> x* = [2, 3]
    #       x_2 >= 3
    q = [1; 1.];
    P = spzeros(2, 2);
    A = spdiagm(0 => ones(2))
    b = [-2.; -3.]
    model = COSMO.Model();
    assemble!(model, P, q, COSMO.Constraint(A, b, COSMO.Nonnegatives), settings = COSMO.Settings(verbose_timing = true, check_termination = 20, ));
    res = COSMO.optimize!(model);

    # update b to give: x_1 >= 0, x_2 >= -1. --> x*=[0, -1]
    @test model.states.KKT_FACTORED
    COSMO.update!(model, b = [0.; 1.])
    res2 = COSMO.optimize!(model);
    @test norm(res2.x - [0; -1.], 2) < 1e-4
    @test res.times.init_factor_time == res2.times.init_factor_time

    # warm-starting of x,y should yield correct s = b - Ax and lead to a fast second solve

    # test full emptying of the model
    COSMO.empty_model!(model)
    @test model.states.IS_OPTIMIZED == false
    @test model.states.IS_ASSEMBLED == false
    @test model.states.KKT_FACTORED == false
    @test model.states.IS_SCALED == false

    assemble!(model, P, q, COSMO.Constraint(A, b, COSMO.Nonnegatives), settings = COSMO.Settings(verbose_timing = true, check_termination = 20));
    @test model.states.IS_ASSEMBLED
    res3 = COSMO.optimize!(model);
    @test abs(res.obj_val - res3.obj_val) < 1e-3
    @test model.states.IS_OPTIMIZED
    @test model.states.KKT_FACTORED
    @test model.states.IS_SCALED
end


# end
