The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "../../../examples/qp.jl"
```

# Quadratic Program

We want to solve the following quadratic program with decision variable `x`:
```math
\begin{array}{ll} \text{minimize} &  1/2 x^\top P x + q^\top x \\
\text{subject to} &  l \leq A x \leq u
\end{array}
```
The problem can be solved with `COSMO` in the following way. Start by defining the problem data

````@example qp
using COSMO, SparseArrays, LinearAlgebra, Test

q = [1; 1.];
P = sparse([4. 1; 1 2]);
A = [1. 1; 1 0; 0 1];
l = [1.; 0; 0];
u = [1; 0.7; 0.7];
nothing #hide
````

First, we decide to solve the problem with two one-sided constraints using `COSMO.Nonnegatives` as the convex set:

````@example qp
Aa = [-A; A]
ba = [u; -l]
constraint1 = COSMO.Constraint(Aa, ba, COSMO.Nonnegatives);
nothing #hide
````

Next, we define the settings object, the model and then assemble everything:

````@example qp
settings = COSMO.Settings(verbose=true);
model = COSMO.Model();
assemble!(model, P, q, constraint1, settings = settings);
res = COSMO.optimize!(model);
nothing #hide
````

Alternatively, we can also use two-sided constraints with `COSMO.Box`:

````@example qp
constraint1 = COSMO.Constraint(A, zeros(3), COSMO.Box(l, u));

model = COSMO.Model();
assemble!(model, P, q, constraint1, settings = settings);
res_box = COSMO.optimize!(model);
nothing #hide
````

Let's check that the solution is correct:

````@example qp
@testset "QP Problem" begin
  @test norm(res.x[1:2] - [0.3; 0.7], Inf) < 1e-3
  @test norm(res_box.x[1:2] - [0.3; 0.7], Inf) < 1e-3
  @test abs(res.obj_val - 1.88) < 1e-3
  @test abs(res_box.obj_val - 1.88) < 1e-3
end
nothing
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

