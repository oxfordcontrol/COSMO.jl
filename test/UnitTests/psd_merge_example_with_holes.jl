
using COSMO, Random, Test, LinearAlgebra, SparseArrays, Random

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

rng = Random.MersenneTwister(2222)

# Generate the problem data
A1 = rand(rng, 9, 9)
A1[1, 5:9] .= 0
A1[2, 6:9] .= 0
A1[3, 6:9] .= 0
A1[4, 9] = 0
A1 = Symmetric(A1, :U)
a1 = vec(A1)

S1 = generate_pos_def_matrix(rng, 9, 0.1, 2)
apply_pattern!(S1, A1)
A = hcat(A1[:])
s = S1[:]
x = rand(rng, 1)
b = A * x + s
B = reshape(b, 9, 9)

Y = generate_pos_def_matrix(rng, 9,  0.1, 1)
y = vec(Y)
P = sparse(zeros(1, 1))
q = -P * x - A' * y

Ct = [COSMO.PsdConeTriangle(45)];
At = zeros(45)
bt = zeros(45)
COSMO.extract_upper_triangle!(A1, At, sqrt(2))
COSMO.extract_upper_triangle!(B, bt, sqrt(2))
At = hcat(At)


settings = COSMO.Settings(decompose = true, eps_abs = 1e-5, eps_rel = 1e-5, max_iter = 5000, complete_dual = true, verbose = true, merge_strategy = COSMO.PairwiseMerge);
model = COSMO.Model();
COSMO.set!(model, P, q, At, bt, Ct, settings);
res = COSMO.optimize!(model);
y = res[1].y
Y = zeros(9, 9)
COSMO.populate_upper_triangle!(Y, y, 1 /sqrt(2))
Y = Symmetric(Y, :U)

Y = sparse(Y)
dropzeros!(Y)
R, C, V = findnz(Y)
C .-= 1
R .-= 1

npzwrite("R.npy", R)
npzwrite("C.npy", C)
npzwrite("V.npy", V)
 #  ws = model
 #  ws.ci = COSMO.ChordalInfo{Float64}(ws.p)

 #  row_ranges = COSMO.get_set_indices(ws.p.C.sets)
 #  C = ws.p.C.sets[1]
 #  psd_row_range = row_ranges[1]
 #  csp = COSMO.find_aggregate_sparsity(ws.p.A, ws.p.b, psd_row_range, C)

 #  ordering = COSMO.find_graph!(ws.ci, csp, C.sqrt_dim, C)
 #  ci = ws.ci
 #  L = ci.L
 #  sntree = COSMO.SuperNodeTree(L, settings.merge_strategy())
 #  #COSMO.merge_cliques!(sntree)

 #  t = sntree
 #  COSMO._merge_cliques!(t, t.strategy)

 #  # since for now we have a graph, not a tree, a post ordering or a parent structure does not make sense. Therefore just number
 #  # the non-empty supernodes in t.snd
 #  t.snd_post = findall(x -> !isempty(x), t.snd)
 #  t.snd_par = -ones(Int64, length(t.snd))

 # #clique_tree_from_graph!(t)
 #  COSMO.clique_intersections!(t.strategy.edges, t.snd)

 #  COSMO.kruskal!(t.strategy.edges, t.num)


 # COSMO.determine_parent_cliques!(t.snd_par, t.snd_child, t.snd, t.post, t.strategy.edges)

 #  t.snd_post = COSMO.post_order(t.snd_par, t.snd_child, t.num)
 #  COSMO.split_cliques!(t.snd, t.sep, t.snd_par, t.snd_post, t.num)


 #  COSMO.reorder_snd_consecutively!(sntree, ordering)





