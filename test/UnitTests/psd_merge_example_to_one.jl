# A unit test with an example problem where you can decompose a 5x5 matrix into two
# 4x4 blocks that overlap so much, that you should merge them back into a 5x5 block.

using COSMO, Random, Test, LinearAlgebra, SparseArrays, Random
rng = Random.MersenneTwister(2222)

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
A1 = rand(rng, 5, 5)
A1 = 0.5 * (A1 + A1')
A1[1, 3] = A1[3, 1] = 0
a1 = vec(A1)

S1 = generate_pos_def_matrix(rng, 5, 0.1, 2)
apply_pattern!(S1, A1)
A = hcat(A1[:])
s = S1[:]
x = rand(rng, 1)
b = A * x + s
B = reshape(b, 5, 5)

Y = generate_pos_def_matrix(rng, 5,  0.1, 1)
y = vec(Y)
P = sparse(zeros(1, 1))
q = -P * x - A' * y

Ct = [COSMO.PsdConeTriangle(15)];
At = zeros(15)
bt = zeros(15)
COSMO.extract_upper_triangle!(A1, At, sqrt(2))
COSMO.extract_upper_triangle!(B, bt, sqrt(2))
At = hcat(At)


settings = COSMO.Settings(decompose = true, complete_dual = false, verbose = true, merge_strategy = COSMO.PairwiseMerge);
model = COSMO.Model();
COSMO.set!(model, P, q, At, bt, Ct, settings);
res = COSMO.optimize!(model);










