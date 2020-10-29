The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "<unknown>/../examples/max_k_cut.jl"
```

Maximum-k-cut problem

```@example max_k_cut
using LinearAlgebra, Random, SparseArrays, COSMO, JuMP

rng = Random.MersenneTwister(1)
n = 4;
W = zeros(n, n);
W[1, 2] = 1; W[1, 4] = 8;
W[2, 3] = 2; W[2, 4] = 10;
W[3, 4] = 6;
W = Symmetric(W)

# Compute the Laplacian matrix of the graph
L = diagm(0 => W * ones(n)) - W;
nothing #hide
```

primal: maximum-cut

```@example max_k_cut
model = JuMP.Model(with_optimizer(COSMO.Optimizer));
@variable(model, Y[1:n, 1:n], PSD);
@objective(model, Max, 1 / 4 * dot(L, Y));
@constraint(model, [i = 1:n], Y[i, i] == 1.);
JuMP.optimize!(model)
```

dual: maximum-cut

```@example max_k_cut
model_dual = JuMP.Model(with_optimizer(COSMO.Optimizer, complete_dual = true));
@variable(model_dual, γ[1:n]);
@objective(model_dual, Min,  sum(γ));
@constraint(model_dual, lmi, Symmetric(-1/4 * L + diagm(γ)) in JuMP.PSDCone());
JuMP.optimize!(model_dual)
obj_val = JuMP.objective_value(model_dual)
Yopt = dual(lmi)
```

primal: maximum-k-cut

```@example max_k_cut
k = 3
model = JuMP.Model(with_optimizer(COSMO.Optimizer));
@variable(model, Y[1:n, 1:n], PSD);
@objective(model, Max, (k - 1) / (2 * k) * dot(L, Y));
@constraint(model, [i = 1:n], Y[i, i] == 1.);
@constraint(model, Y[1, 2] >= -1/(k-1));
@constraint(model, Y[1, 4] >= -1/(k-1));
@constraint(model, Y[2, 3] >= -1/(k-1));
@constraint(model, Y[2, 4] >= -1/(k-1));
@constraint(model, Y[3, 4] >= -1/(k-1));
JuMP.optimize!(model)
Yopt = JuMP.value.(Y);
nothing #hide
```

dual: maximum-k-cut

```@example max_k_cut
model_dual = JuMP.Model(with_optimizer(COSMO.Optimizer, decompose = false))#complete_dual = true));
@variable(model_dual, γ[1:n]);
@variable(model_dual, λ[1:6]);
@objective(model_dual, Min,  sum(γ) + 1 / (k-1) * sum(λ) );
@constraint(model_dual, λ .>= 0);
Δ = [0. λ[1] 0 λ[2]; λ[1] 0 λ[3] λ[4]; 0 λ[3] 0 λ[5]; λ[2] λ[4] λ[5] 0 ]
Δ = [0. λ[1] λ[2] λ[3]; 0 0 λ[4] λ[5]; 0 0 0 λ[6]; 0 0 0 0 ]
@constraint(model_dual, lmi, Symmetric(-(k - 1) / (2 * k) * L + diagm(γ) - Δ) in JuMP.PSDCone());
JuMP.optimize!(model_dual)
obj_val = JuMP.objective_value(model_dual)
```

```@example max_k_cut
Yopt = JuMP.value.(Y)
γopt = JuMP.value.(γ)
λopt = JuMP.value.(λ)
obj_val = JuMP.objective_value(model_dual)



factor = cholesky(Yopt, Val(true); check = false);
V = Matrix((factor.P * factor.L)');
# normalize columns
for i in 1:n
    V[:, i] ./= norm(V[:, i]);
end
```

Rounding after Frieze and Jerrum
choose k random vectors gj from N(0, 1)

```@example max_k_cut
G = randn(rng, n, k)
y = zeros(Int64, n)
for i in 1:n
    val = -Inf
    for j = 1:k
        prod = dot(G[:, j], V[:, i])
        if prod >= val
            val = prod
            y[i] = j
        end
    end
end
```

or each vertexi∈V, consider thekdot products of vectorviwitheach of thekrandom vectors,vi·g1,vi·g2,...,vi·gk. One of these dot products is # maximum.Assign the vertex the label of the random vector with which it has the maximum dot product

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

