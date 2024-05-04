# # Arbitrary Precision
#
# COSMO allows you to solve problems with arbitrary floating-point precision, e.g. by using `BigFloat` problem data. To do this, the desired floating point type has to be consistent across the model `COSMO.Model{<: AbstractFloat}`, the input data and the (optional) settings object `COSMO.Settings{<: AbstractFloat}`.
# As an example, assume we want to solve the following quadratic program:
#
# $$
# \begin{array}{ll} \text{minimize} &  1/2 x^\top P x + q^\top x \\
# \text{subject to} &  l \leq A x \leq u
# \end{array}
# $$
# where $P = \begin{bmatrix} 4 & 1 \\ 1 & 2\end{bmatrix}$, $q = [1, 1]^\top$, $A = \begin{bmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1\end{bmatrix}$ and $l= [1,0, 0]^\top$, $u=[1, 0.7, 0.7]^\top$. We start by creating the model with the desired precision:
using COSMO, LinearAlgebra, SparseArrays
model = COSMO.Model{BigFloat}()

#-
# Next, we define the problem data as `BigFloat` arrays and create the constraint:
q = BigFloat[1; 1.]
P = sparse(BigFloat[4. 1; 1 2])
A = BigFloat[1. 1; 1 0; 0 1]
l = BigFloat[1.; 0; 0]
u = BigFloat[1; 0.7; 0.7]

constraint = COSMO.Constraint(A, zeros(BigFloat, 3), COSMO.Box(l, u))

# Notice that the constraint type parameter is dependent on the input data. The same is true for the constraint set `Box`. Next, we define the settings
settings = COSMO.Settings{BigFloat}(verbose = true, kkt_solver = QdldlKKTSolver)

# and assemble and solve the problem:
assemble!(model, P, q, constraint, settings = settings)
result = COSMO.optimize!(model);


# Moreover, notice that when no type parameter is specified, all objects default to `Float64`:
model = COSMO.Model()

# Two limitations to arbitrary precision:
# - Since we call the LAPACK function `syevr` for eigenvalue decompositions, we currently only support solving problems with PSD constraints in `Float32` and `Float64`.
# - We suggest to use the pure Julia QDLDL linear system solver (`kkt_solver = QdldlKKTSolver`) when working with arbitrary precision types as some of the other available solvers don't support all available precisions.

#md # !!! note
#md #    If you want to use `COSMO` directly with `MathOptInterface`, you can use: `COSMO.Optimizer{<: AbstractFloat}` as your optimizer. Again, the problem data precision of your MathOptInterface-model has to agree with the optimizer's precision.
