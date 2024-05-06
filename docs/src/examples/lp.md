The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "../../../examples/lp.jl"
```

# Linear Program

We want to solve the following linear program with decision variable `x`:
```math
\begin{array}{ll} \text{minimize} &  c^\top x\\
\text{subject to} &  A x \leq b \\
                  &  x \geq 1 \\
                  &  x_2 \geq 5 \\
                  &  x_1 + x_3 \geq 4.
\end{array}
```
The problem can be solved with `COSMO` in the following way:

````@example lp
using COSMO, LinearAlgebra, SparseArrays, Test
````

````@example lp
##Define problem data:
c = [1; 2; 3; 4.];
A = Matrix(1.0I, 4, 4);
b = [10.; 10; 10; 10];
n = 4;
nothing #hide
````

Create the constraints $Ax + b \in \mathcal{K}$:

````@example lp
# Ax <= b
c1 = COSMO.Constraint(-A, b, COSMO.Nonnegatives);
# x >= 1
c2 = COSMO.Constraint(Matrix(1.0I, n, n), -ones(n), COSMO.Nonnegatives);
# x2 >= 5
c3 = COSMO.Constraint(1, -5, COSMO.Nonnegatives, n, 2:2);
# x1 + x3 >= 4
c4 = COSMO.Constraint([1 0 1 0], -4, COSMO.Nonnegatives);
nothing #hide
````

Define matrix $P$ and vector $q$ for the objective function:

````@example lp
P = spzeros(4, 4);
q = c;
nothing #hide
````

Create, assemble and solve the model:

````@example lp
settings = COSMO.Settings(verbose=true, eps_abs = 1e-4, eps_rel = 1e-5);
model = COSMO.Model();
assemble!(model, P, q, [c1; c2; c3; c4], settings = settings);
res = COSMO.optimize!(model)
````

Compare the result to the known solution:

````@example lp
@test isapprox(res.x[1:4], [3; 5; 1; 1], atol=1e-2, norm = (x -> norm(x, Inf)))
````

````@example lp
@test isapprox(res.obj_val, 20.0, atol=1e-2)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

