The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "../../../examples/closest_correlation_matrix.jl"
```

# Closest Correlation Matrix

We consider the problem of finding the closest correlation matrix $X$ to a given random matrix $C$.
With closest correlation matrix we mean a positive semidefinite matrix with ones on the diagonal.
The problem is given by:
```math
\begin{array}{ll} \text{minimize} &  \frac{1}{2}||X - C||_F^2\\
\text{subject to} &  X_{ii} = 1, \quad i=1,\dots,n \\
                  &  X \succeq 0.
\end{array}
```
Notice that we use `JuMP` to model the problem. `COSMO` is chosen as the backend solver. And COSMO-specific settings are passed using
the `optimizer_with_attributes()` function.

````@example closest_correlation_matrix
using COSMO, JuMP, LinearAlgebra, SparseArrays, Test, Random
````

````@example closest_correlation_matrix
rng = Random.MersenneTwister(12345);
# create a random test matrix C
n = 8;
C = -1 .+ rand(rng, n, n) .* 2;
c = vec(C);
nothing #hide
````

Define problem in `JuMP`:

````@example closest_correlation_matrix
q = -vec(C);
r = 0.5 * vec(C)' * vec(C);
m = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "verbose" => true));
@variable(m, X[1:n, 1:n], PSD);
x = vec(X);
@objective(m, Min, 0.5 * x' * x  + q' * x + r);
for i = 1:n
  @constraint(m, X[i, i] == 1.);
end
````

Solve the `JuMP` model with `COSMO` and query the solution `X_sol`:

````@example closest_correlation_matrix
status = JuMP.optimize!(m);
obj_val = JuMP.objective_value(m);
X_sol = JuMP.value.(X);
````

Double check result against known solution:

````@example closest_correlation_matrix
known_opt_val = 8.75427
known_solution =  [1.0         0.354428    0.238389    0.323015  -0.117854    0.620975    0.0969092   0.231657
  0.354428    1.0        -0.590998   -0.172179  -0.571414    0.232656    0.109805    0.0229169
  0.238389   -0.590998    1.0         0.379909   0.0205233  -0.237496    0.0431976  -0.0364916
  0.323015   -0.172179    0.379909    1.0       -0.250821    0.349762    0.13669     0.49255
 -0.117854   -0.571414    0.0205233  -0.250821   1.0         0.145922   -0.0418897   0.0604382
  0.620975    0.232656   -0.237496    0.349762   0.145922    1.0         0.457861    0.0215352
  0.0969092   0.109805    0.0431976   0.13669   -0.0418897   0.457861    1.0        -0.626215
  0.231657    0.0229169  -0.0364916   0.49255    0.0604382   0.0215352  -0.626215    1.0   ];

@test isapprox(obj_val, known_opt_val , atol=1e-3)
````

````@example closest_correlation_matrix
@test norm(X_sol - known_solution, Inf) < 1e-3
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

