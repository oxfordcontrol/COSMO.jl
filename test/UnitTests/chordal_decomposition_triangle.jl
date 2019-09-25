using COSMO, Random, Test, LinearAlgebra, SparseArrays, Random
rng = Random.MersenneTwister(144545)


# We create a test problem with 4 convex sets (2 of them decomposable), 1 zero set in the middle

# define convex set
C = [COSMO.PsdCone(16); COSMO.ZeroSet(2); COSMO.PsdCone(16); COSMO.PsdCone(9)]

m = 1

# CONSTRAINT 1
A1 = rand(rng, 4, 4)
A1 = 0.5 * (A1 + A1')
# eliminate certain entries to make matrix chordal
A1[1, 3] = A1[1, 4] = A1[3, 1] = A1[4, 1] = 0
a1 = vec(A1)
# find solution matrix S that is pos def
S1 = A1 + (eigvals(A1)[1] + 1) * Matrix(1.0I, 4, 4)

# CONSTRAINT 2
a2 = rand(rng, 2, 1)
s2 = [0; 0]

# CONSTRAINT 3
A3 = rand(rng, 4, 4)
A3 = 0.5 * (A3 + A3')
# eliminate certain entries to make matrix chordal
A3[2, 4] = A3[4, 2] =  A3[1, 3]  = A3[3, 1] = A3[2, 3] = A3[3, 2] = 0
a3 = vec(A3)
S3 = A3 + (eigvals(A3)[1] + 1) * Matrix(1.0I, 4, 4)

# CONTRAINT 4 (dense and undecomposable)
A4 = rand(rng, 3, 3)
A4 = 0.5 * (A4 + A4')
a4 = vec(A4)
S4 = A4 + (eigvals(A4)[1] + 1) * Matrix(1.0I, 3, 3)

# Find b such that problem is feasible
x = rand(rng, 1)
A = [a1; a2; a3; a4]
s = [vec(S1); s2; vec(S3); vec(S4)]
b = A * x + s

# Find q such that problem is feasible
Y1 = generate_pos_def_matrix(rng, 4,  0.1, 1)
y2 = rand(rng, 2)
Y3 = generate_pos_def_matrix(rng, 4,  0.1, 1)
Y4 = generate_pos_def_matrix(rng, 3,  0.1, 1)
y = [vec(Y1); y2; vec(Y3); vec(Y4)]

P = sparse(zeros(1, 1))
q = -P * x - A' * y

# --------------------------
# CONFIG 1: PsdSquare + no decomposition
# --------------------------

model = COSMO.Model()
settings = COSMO.Settings(decompose = false)
COSMO.set!(model, P, q, A, b, C, settings)
res1 = COSMO.optimize!(model);
# Y1_sol1 = reshape(res1.y[1:16], 4, 4)
# Y2_sol1 = reshape(res1.y[19:34], 4, 4)
S1_sol1 = reshape(res1.s[1:16], 4, 4)
S2_sol1 = reshape(res1.s[19:34], 4, 4)


# --------------------------
# CONFIG 2: PsdSquare +  decomposition
# --------------------------

model = COSMO.Model()
settings = COSMO.Settings(decompose = true)
COSMO.set!(model, P, q, A, b, C, settings)
res2 = COSMO.optimize!(model);
S1_sol2 = reshape(res2.s[1:16], 4, 4)
S2_sol2 = reshape(res2.s[19:34], 4, 4)


# --------------------------
# Only keep upper triangular data of the problem
# --------------------------

# Notice that we have to rescale the off-diagonal entries of s by sqrt(2) since our projection assumes
# that
Ct = [COSMO.PsdConeTriangle(10); COSMO.ZeroSet(2); COSMO.PsdConeTriangle(10); COSMO.PsdConeTriangle(6)]

a1t = zeros(10)
COSMO.extract_upper_triangle!(A1, a1t, sqrt(2))
a3t = zeros(10)
COSMO.extract_upper_triangle!(A3, a3t, sqrt(2))
a4t = zeros(6)
COSMO.extract_upper_triangle!(A4, a4t, sqrt(2))
At = [a1t; a2; a3t; a4t]

B1 = reshape(b[1:16], 4, 4)
B3 = reshape(b[19:34], 4, 4)
B4 = reshape(b[35:end], 3, 3)
b1t = zeros(10)
b3t = zeros(10)
b4t = zeros(6)
COSMO.extract_upper_triangle!(B1, b1t, sqrt(2))
COSMO.extract_upper_triangle!(B3, b3t, sqrt(2))
COSMO.extract_upper_triangle!(B4, b4t, sqrt(2))
bt = [b1t; b[17:18]; b3t; b4t]


# --------------------------
# CONFIG 3: PsdTriangle + no decomposition
# --------------------------

model = COSMO.Model()
settings = COSMO.Settings(decompose = false)
COSMO.set!(model, P, q, At, bt, Ct, settings)
res3 = COSMO.optimize!(model);
S1_sol3 = zeros(4, 4)
COSMO.populate_upper_triangle!(S1_sol3, res3.s[1:10], 1. / sqrt(2))
S1_sol3 = Symmetric(S1_sol3)
S2_sol3 = zeros(4, 4)
COSMO.populate_upper_triangle!(S2_sol3, res3.s[13:22], 1. / sqrt(2))
S2_sol3 = Symmetric(S2_sol3)

# --------------------------
# CONFIG 4: PsdTriangle +  decomposition
# --------------------------

model = COSMO.Model()
settings = COSMO.Settings(decompose = true)
COSMO.set!(model, P, q, At, bt, Ct, settings)
res4 = COSMO.optimize!(model);
S1_sol4 = zeros(4, 4)
COSMO.populate_upper_triangle!(S1_sol4, res4.s[1:10], 1. / sqrt(2))
S1_sol4 = Symmetric(S1_sol4)
S2_sol4 = zeros(4, 4)
COSMO.populate_upper_triangle!(S2_sol4, res4.s[13:22], 1. / sqrt(2))
S2_sol4 = Symmetric(S2_sol4)

@testset "Decomposition with PSDTriangle" begin
  @test abs(res1.obj_val - res2.obj_val) < 1e-4
  @test abs(res1.obj_val - res3.obj_val) < 1e-4
  @test abs(res2.obj_val - res4.obj_val) < 1e-4
  @test norm(S1_sol1 - S1_sol2, Inf) < 1e-8
  @test norm(S1_sol2 - S1_sol3, Inf) < 1e-8
  @test norm(S1_sol3 - S1_sol4, Inf) < 1e-8
  @test norm(S2_sol1 - S2_sol2, Inf) < 1e-8
  @test norm(S2_sol2 - S2_sol3, Inf) < 1e-8
  @test norm(S2_sol3 - S2_sol4, Inf) < 1e-8
end
nothing


