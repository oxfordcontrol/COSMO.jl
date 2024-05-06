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
nothing #hide
````

Double check result against known solution:

````@example closest_correlation_matrix
known_opt_val = 12.5406
known_solution =  [
  1.0         0.732562   -0.319491   -0.359985   -0.287543   -0.15578     0.0264044  -0.271438;
  0.732562    1.0         0.0913246  -0.0386357   0.299199   -0.122733    0.126612   -0.187489;
 -0.319491    0.0913246   1.0        -0.0863377   0.432948    0.461783   -0.248641   -0.395299;
 -0.359985   -0.0386357  -0.0863377   1.0         0.503379    0.250601    0.141151    0.286088;
 -0.287543    0.299199    0.432948    0.503379    1.0        -0.0875199   0.137518    0.0262425;
 -0.15578    -0.122733    0.461783    0.250601   -0.0875199   1.0        -0.731556    0.0841783;
  0.0264044   0.126612   -0.248641    0.141151    0.137518   -0.731556    1.0        -0.436274;
 -0.271438   -0.187489   -0.395299    0.286088    0.0262425   0.0841783  -0.436274    1.0  ];
@test isapprox(obj_val, known_opt_val , atol=1e-3)
````

````@example closest_correlation_matrix
@test norm(X_sol - known_solution, Inf) < 1e-3
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

