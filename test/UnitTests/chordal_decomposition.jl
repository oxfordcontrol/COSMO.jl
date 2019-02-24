using COSMO, SparseArrays, LinearAlgebra, FileIO, Test

if !@isdefined test_problem_path
  test_problem_path = "./testproblems/";
end

data = load(test_problem_path * "decomp_testProblem.jld2")

P = data["P"];
q = data["q"];
A = data["A"];
b = data["b"];
obj_true = data["objTrue"];
cs = COSMO.Constraint(-A, b, COSMO.PsdCone)

model = COSMO.Model()
settings = COSMO.Settings(decompose = false, obj_true = obj_true)
assemble!(model, P, q, [cs], settings)
# # solve problem with chordal decomposition turned off
res = COSMO.optimize!(model);


P = data["P"];
q = data["q"];
A = data["A"];
b = data["b"];
cs = COSMO.Constraint(-A, b, COSMO.PsdCone)
# # solve problem with sparsity exploitation
model = COSMO.Model()
settings_decomp = COSMO.Settings(decompose = true, obj_true = obj_true, verbose_timing = true)
assemble!(model, P, q, [cs], settings_decomp);
res_decomp = COSMO.optimize!(model);


#  compare results
@testset "Chordal Decomposition" begin
  @test abs(res_decomp.obj_val - res.obj_val) < 1e-3
  @test abs(res_decomp.obj_val - obj_true) < 1e-3
  @test norm(res_decomp.x - res.x,Inf) < 1e-3
end
nothing