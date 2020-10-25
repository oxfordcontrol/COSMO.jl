# test some utility functions specific to the python interface cosmo-python
using COSMO, Test, LinearAlgebra, Random, SparseArrays

rng = Random.MersenneTwister(41)


@testset "cosmo-python related utility function" begin

    # create a set of convex sets from a dictionary
    cone = Dict("f" => 2, "l" => 3, "q" => [3; 4], "s" => [3; 6; 10], "ep" => 2, "ed" => 1, "p" => [0.3; -0.4])
    convex_sets = COSMO.convex_sets_from_dict(cone)
    @test convex_sets[3] isa COSMO.SecondOrderCone
    @test convex_sets[3].dim == 3
    @test convex_sets[6] isa COSMO.PsdConeTriangle
    @test convex_sets[11] isa COSMO.PowerCone
    @test convex_sets[11].α == 0.3


    # solve QP using set! function with sparse matrix input in (rowval, colptr, nzval)-format
    q = [1; 1.];
    P = sparse([4. 1; 1 2]);
    A1 = [1. 1; 1 0; 0 1];
    A = sparse([A1; -A1])

    l = [1.; 0; 0];
    u = [1; 0.7; 0.7];
    b = [u; -l]
    cone = Dict("l" => 6)
    m, n = size(A)
    model = COSMO.Model();
    COSMO.set!(model, P.rowval, P.colptr, P.nzval, q, A.rowval, A.colptr, A.nzval, b, cone, m, n);
    res = COSMO.optimize!(model);
end
