The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "<unknown>/../examples/sum_abs_k_eigenvalues.jl"
```

# Minimizing the sum of the k-largest λ

We show how to find the sum of absolute value of the `k` largest eigenvalues of a symmetric matrix $A \in \mathbb{S}^n$. This problem can be solved as a semidefinite program. The primal and dual forms are stated in *Alizadeh* \[1\]:
```math
\begin{array}{llll} \text{maximize} &   \text{Tr}(AY) - \text{Tr}(AW) &    \text{minimize} &  kz + Tr(U) + Tr(V)            \\
\text{subject to} &  \text{Tr}(Y + W) = k                             &   \text{subject to} &   zI + V - A \succeq 0            \\
                  &  0 \preceq  Y \preceq I                           &  &                    zI + U + A \succeq 0             \\
                  &  0 \preceq  W \preceq I                           & &                     U, V \succeq 0,
\end{array}
```
where $Y, W$ are the variables of the primal and $U, V$ are the variables of the dual problem.

```@example sum_abs_k_eigenvalues
using LinearAlgebra, JuMP, COSMO, Random
rng = Random.MersenneTwister(212)

n = 10
A = 5 .* randn(rng, 10, 10)
A = Symmetric(A, :U)
```

We are interested in minimizing the sum of absolute values of the `k=3` largest eigenvalues. Let's formulate the problem in `JuMP` with `COSMO` as the backend solver:

```@example sum_abs_k_eigenvalues
k = 3
model = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose=true));
@variable(model, Y[1:n, 1:n], PSD);
@variable(model, W[1:n, 1:n], PSD);

@objective(model, Max, tr(A * Y) - tr(A * W));
@constraint(model, tr(Y + W) == k);
@constraint(model, Symmetric(I - Y) in PSDCone());
@constraint(model, Symmetric(I - W) in PSDCone());
status = JuMP.optimize!(model)
```

```@example sum_abs_k_eigenvalues
opt_objective = JuMP.objective_value(model)
```

Now, we can check the solution by computing the sum of the absolute value of the 3-largest eigenvalues:

```@example sum_abs_k_eigenvalues
k_λ_abs = sum(sort(abs.(eigen(A).values), rev = true)[1:k])
```

### Solve the dual

Alternatively, we can solve the dual problem:

```@example sum_abs_k_eigenvalues
model = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose=true));
@variable(model, V[1:n, 1:n], PSD);
@variable(model, U[1:n, 1:n], PSD);
@variable(model, z);

@objective(model, Min, k * z + tr(V) + tr(U));
@constraint(model, Symmetric(z .* diagm(0 => ones(n)) + V - A) in PSDCone());
@constraint(model, Symmetric(z .* diagm(0 => ones(n)) + U + A) in PSDCone());
status = JuMP.optimize!(model)
```

```@example sum_abs_k_eigenvalues
opt_objective = JuMP.objective_value(model)
```

This gives the same result.
## Problem with A as variable

Above problems are mostly helpful for illustrative purpose. It is obviously easier to find the sum of the k-largest eigenvalues by simply computing the eigenvalues of $A$. However, above results become useful if finding $A$ itself is part of the problem. For example, assume we want to find a valid matrix $A$ under the constraints: $C\, \text{vec}(A) = b$ with the minimum sum of absolute values of the k-largest eigenvalues. We can then solve the equivalent problem:
```math
\begin{array}{ll} \text{minimize} &  kz + Tr(U) + Tr(V)     \\
 \text{subject to} &   C \text{vec}(A) = b \\
                   & zI + V - A \succeq 0            \\
                   &    zI + U + A \succeq 0             \\
                   &   U, V \succeq 0.
\end{array}
```

## References
[1] Alizadeh - Interior point methods in semidefinite programming with applications to combinatorial optimization (1995)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

