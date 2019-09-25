# This unit test checks the positive semidefinite completion routine of the dual variable
# We create a random feasible problem with one 4x4 PSD constraint that has a chordal structure
# and can be decomposed into two smaller PSD constraints
# We consider 2 cases:
# 1. Solve SDP  without decomposing the constraint (both primal slack variable S and dual variable Y should be positive semidefinite)
# 2.    -"-     with          -"-                 and with  completing the dual variable Y  (Y should be pos semdef)

using COSMO, Random, Test, LinearAlgebra, SparseArrays, Random
rng = Random.MersenneTwister(144545)


# Generate the problem data
C = [COSMO.PsdCone(16)]
m = 1

A1 = rand(rng, 4, 4)
A1 = 0.5 * (A1 + A1')
A1[1, 3] = A1[1, 4] = A1[3, 1] = A1[4, 1] = 0
a1 = vec(A1)

S1 = generate_pos_def_matrix(rng, 4, 0.1, 2)
apply_pattern!(S1, A1)
A = hcat(A1[:])
s = S1[:]
x = rand(rng, 1)
b = A * x + s

Y = generate_pos_def_matrix(rng, 4,  0.1, 1)
y = vec(Y)
P = sparse(zeros(1, 1))
q = -P * x - A' * y


# Run two test cases
cases = [
  COSMO.Settings(decompose = false, verbose_timing = true);
  COSMO.Settings(decompose = true, complete_dual = true, verbose_timing = true, merge_strategy = COSMO.NoMerge)
  ]
results = Array{COSMO.Result{Float64}}(undef, 2);

for i = 1:2
  model = COSMO.Model()
  settings = cases[i]
  COSMO.set!(model, P, q, A, b, C, settings)
  res = COSMO.optimize!(model);
  results[i] = res
end

@testset "PSD Completion (PsdConeSquare)" begin
  Y_sol1 = reshape(results[1].y, 4, 4)
  Y_sol1 = 0.5 * (Y_sol1 + Y_sol1')
  Y_sol3 = reshape(results[2].y, 4, 4)
  @test minimum(eigvals(Y_sol1)) > -1e-6
  @test minimum(eigvals(Y_sol3)) > -1e-6

end

# PSDCone Triangle
# We consider the same cases but now work only with the upper triangular data. The results should be the same
Ct = [COSMO.PsdConeTriangle(10)];
B = reshape(b, 4, 4)
At = zeros(10)
bt = zeros(10)
COSMO.extract_upper_triangle!(A1, At, sqrt(2))
COSMO.extract_upper_triangle!(B, bt, sqrt(2))
At = hcat(At)
results_triangle = Array{COSMO.Result{Float64}}(undef, 3);


cases_triangle = [
  COSMO.Settings(decompose = false);
  COSMO.Settings(decompose = true, complete_dual = true, merge_strategy = COSMO.NoMerge)
  ]
for i = 1:2
  model_tri = COSMO.Model()
  settings = cases_triangle[i]
  COSMO.set!(model_tri, P, q, At, bt, Ct, settings)
  res = COSMO.optimize!(model_tri);
  results_triangle[i] = res
end

@testset "PSD Completion (PsdConeTriangle)" begin
  Y_sol1t = matrixify(results_triangle[1].y)
  Y_sol2t = matrixify(results_triangle[2].y)
  @test minimum(eigvals(Y_sol1t)) > -1e-6
  @test minimum(eigvals(Y_sol2t)) > -1e-6
end







