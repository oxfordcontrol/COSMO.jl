# Examples
Some example problems that are solved with COSMO can be found in the [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples) folder.

In the first example the native solver interface is used to define and solve a Linear Program. In the second example `JuMP` is used to describe the problem and COSMO is set as the solver backend.

## Linear Program
We want to solve the following linear program with decision variable `x`:
```math
\begin{array}{ll} \mbox{minimize} &  c^\top x\\
\mbox{subject to} &  A x \leq b \\
                  &  x \geq 1 \\
                  &  x_2 \geq 5 \\
                  &  x_1 + x_3 \geq 4.
\end{array}
```
The problem can be solved with COSMO in the following way:

```julia
using COSMO, LinearAlgebra, SparseArrays, Test

c = [1; 2; 3; 4.]
A = Matrix(1.0I, 4, 4)
b = [10; 10; 10; 10]
n = 4
# -------------------
# create constraints A * x + b in set
# -------------------
# Ax <= b
c1 = COSMO.Constraint(-A, b, COSMO.Nonnegatives)
# x >= 1
c2 = COSMO.Constraint(Matrix(1.0I, n, n), -ones(n), COSMO.Nonnegatives)
# x2 >= 5
c3 = COSMO.Constraint(1, -5, COSMO.Nonnegatives, n, 2:2)
# x1 + x3 >= 4
c4 = COSMO.Constraint([1 0 1 0], -4, COSMO.Nonnegatives)

# -------------------
# define cost function
# -------------------
P = spzeros(4, 4)
q = c

# -------------------
# assemble solver model
# -------------------
settings = COSMO.Settings(max_iter=2500, verbose=true, eps_abs = 1e-4, eps_rel = 1e-5)
model = COSMO.Model()
assemble!(model, P, q, [c1; c2; c3; c4], settings = settings)
res = COSMO.optimize!(model);

@testset "Linear Problem" begin
  @test isapprox(res.x[1:4], [3; 5; 1; 1], atol=1e-2, norm = (x -> norm(x, Inf)))
  @test isapprox(res.obj_val, 20.0, atol=1e-2)
end
```

## Closest Correlation Matrix
We consider the problem of finding the closest correlation matrix `X` to a given random matrix `C`. With closest correlation matrix we mean a positive semidefinite matrix with ones on the diagonal. The problem is given by:
```math
\begin{array}{ll} \mbox{minimize} &  \frac{1}{2}||X - C||_F^2\\
\mbox{subject to} &  X_{ii} = 1, \quad i=1,\dots,n \\
                  &  X \succeq 0.
\end{array}
```
Notice how `JuMP` is used to describe the problem. COSMO is chosen as the backend solver using JuMP's `with_optimizer()` function.
```julia
using COSMO, JuMP, LinearAlgebra, SparseArrays, Test, Random
rng = Random.MersenneTwister(12345);

# create a random test matrix C
n = 8
C = -1 .+ rand(rng, n, n) .* 2;
c = vec(C);

# define problem in JuMP
q = -vec(C);
r = 0.5 * vec(C)' * vec(C);
m = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose=true, eps_abs = 1e-4));
@variable(m, X[1:n, 1:n], PSD);
x = vec(X);
@objective(m, Min, 0.5 * x' * x  + q' * x + r)
for i = 1:n
  @constraint(m, X[i, i] == 1.)
end

# solve and get results
status = JuMP.optimize!(m)
obj_val = JuMP.objective_value(m)
X_sol = JuMP.value.(X)
```

## Logistic Regression
Logistic regression problems can be solved using exponential cone constraints. An example on how to use COSMO to solve a logistic regression problem is presented in [/examples/logistic\_regression\_regularization.ipynb](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/logistic_regression_regularization.ipynb).
