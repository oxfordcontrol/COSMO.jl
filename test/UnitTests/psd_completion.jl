# This unit test checks the positive semidefinite completion routine of the dual variable
# We create a random feasible problem with one 4x4 PSD constraint that has a chordal structure
# and can be decomposed into two smaller PSD constraints
# We consider 3 cases:
# 1. Solve SDP  without decomposing the constraint (both primal slack variable S and dual variable Y should be positive semidefinite)
# 2.    -"-     with    decomposition but without completing the dual variable Y (Y will not necessarily be pos semdef)
# 3.    -"-     -"-          -"-      and with  completing the dual variable Y  (Y should be pos semdef)

using COSMO, Random, Test, LinearAlgebra, SparseArrays, Random
rng = Random.MersenneTwister(144545)

function generate_pos_def_matrix(rng, n, aMin, aMax)
  X = rand(rng, n, n)
  # any real square matrix can be QP decomposed into a orthogonal matrix and an uppertriangular matrix R
  Q, R = qr(X)
  eigs = rand(rng, n).*(aMax .- aMin) .+ aMin
  X = Q * Diagonal(eigs) * Q'
  X = 0.5 * (X + X')
  return X
end

function apply_pattern!(A, pattern)
  m, n = size(A)
  for i = 1:m, j = 1:n
    if pattern[i, j] == 0
      A[i, j] = 0.
    end
  end
end

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


# Run three test cases
cases = [
  COSMO.Settings(decompose = false, verbose_timing = true);
  COSMO.Settings(decompose = true, complete_dual = false, verbose_timing = true, merge_strategy = COSMO.NoMerge);
  COSMO.Settings(decompose = true, complete_dual = true, verbose_timing = true, merge_strategy = COSMO.NoMerge)
  ]
results = Array{COSMO.Result{Float64}}(undef, 3);

for i = 1:3
  model = COSMO.Model()
  settings = cases[i]
  COSMO.set!(model, P, q, A, b, C, settings)
  res = COSMO.optimize!(model);
  results[i] = res
end

@testset "PSD Completion (PsdConeSquare)" begin
  Y_sol1 = reshape(results[1].y, 4, 4)
  Y_sol1 = 0.5 * (Y_sol1 + Y_sol1')
  Y_sol2 = reshape(results[2].y, 4, 4)
  Y_sol3 = reshape(results[3].y, 4, 4)
  @test minimum(eigvals(Y_sol1)) > -1e-6
  @test minimum(eigvals(Y_sol2)) < -1e-6
  @test minimum(eigvals(Y_sol3)) > -1e-6
  @test norm(Y_sol1 - Y_sol3) < 1e-3

end

# PSDCone Triangle
# We consider the same three cases but now work only with the upper triangular data. The results should be the same
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
  COSMO.Settings(decompose = true, complete_dual = false, merge_strategy = COSMO.NoMerge);
  COSMO.Settings(decompose = true, complete_dual = true, merge_strategy = COSMO.NoMerge)
  ]
for i = 1:3
  model_tri = COSMO.Model()
  settings = cases_triangle[i]
  COSMO.set!(model_tri, P, q, At, bt, Ct, settings)
  res = COSMO.optimize!(model_tri);
  results_triangle[i] = res
end

@testset "PSD Completion (PsdConeTriangle)" begin
  Y_sol1t = zeros(4, 4)
  COSMO.populate_upper_triangle!(Y_sol1t, results_triangle[1].y, 1/sqrt(2))
  Y_sol1t = Symmetric(Y_sol1t)
  Y_sol2t = zeros(4, 4)
  COSMO.populate_upper_triangle!(Y_sol2t, results_triangle[2].y, 1/sqrt(2))
  Y_sol2t = Symmetric(Y_sol2t)
  Y_sol3t = zeros(4, 4)
  COSMO.populate_upper_triangle!(Y_sol3t, results_triangle[3].y, 1/sqrt(2))
  Y_sol3 = Symmetric(Y_sol3t)
  @test minimum(eigvals(Y_sol1t)) > -1e-6
  @test minimum(eigvals(Y_sol2t)) < -1e-6
  @test minimum(eigvals(Y_sol3t)) > -1e-6
  @test abs(results[1].obj_val - results_triangle[1].obj_val) < 1e-3
end







