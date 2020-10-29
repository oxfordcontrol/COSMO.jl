using COSMO, Test, LinearAlgebra, SparseArrays
# this is the problem
# "sdp_sigma_max_atom"
# from the Convex.jl problem set that caused errors in the
# compact_transformation

# max   σ_max(Y)
# s.t.  Y[2, 1] <= 4, Y[2, 2] >= 3, sum(Y) >= 12
# where Y in R_{3x3}
#
# which is equivalent to
#   min t
#   s.t. [t * I Y; Y' t * I]  ⪰ 0

# define x := [t; vec(Y)]
@testset "Chordal LMI-SDP - Maximum singular value problem" begin
    q = [1.0; zeros(9)]
    con1 = COSMO.Constraint([0 0 -1.0 zeros(1, 7)], 4.0, COSMO.Nonnegatives)
    con2 = COSMO.Constraint([0 0 0 0 0 1.0 zeros(1, 4)], -3.0, COSMO.Nonnegatives)
    con3 = COSMO.Constraint([0.0 ones(1, 9)], -12.0, COSMO.Nonnegatives)

    A_lmi1 = [-1.;0; -1; 0; 0; -1; 0; 0; 0; -1; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; -1]
    A_lmi2 = spzeros(21, 9)
    A_lmi2[7, 1] = A_lmi2[8, 2] = A_lmi2[9, 3] = A_lmi2[11, 4] = A_lmi2[12, 5] = A_lmi2[13, 6] = A_lmi2[16, 7] = A_lmi2[17, 8] = A_lmi2[18, 9] = -sqrt(2)
    A_lmi = [A_lmi1 A_lmi2]
    con4 = COSMO.Constraint(-A_lmi, zeros(21), COSMO.PsdConeTriangle(21))

    model = COSMO.Model()
    assemble!(model, spzeros(10, 10), q, [con1; con2; con3; con4], settings = COSMO.Settings(compact_transformation = true, decompose = true, eps_abs = 1e-5, eps_rel = 1e-5))
    res = COSMO.optimize!(model);
    Y = reshape(res.x[2:end], 3, 3)
    # check inequality constraints
    @test Y[2, 1] <= 4
    @test Y[2, 2] >= 3
    @test sum(Y) - 12.0 >= -1e-3
    # check maximum singular value
    σ_max = maximum(svd(Y).S)
    t = res.x[1]
    @test abs(σ_max - t ) <= 1e-3

end
