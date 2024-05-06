The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "../../../examples/svm_primal.jl"
```

# Support Vector Machine
We are showing how to solve a support vector machine problem with COSMO (and JuMP).

## Generating the Dataset
We want to classify the points in this example dataset with $m = 100$ samples and $n = 2$ features:

````@example svm_primal
using Distributions: MvNormal
using Plots, LinearAlgebra, SparseArrays, Random, Test
````

````@example svm_primal
# Generate dataset
rng = Random.MersenneTwister(123);
num_samples = 100;
Xpos = rand(rng, MvNormal([1.5, 1.5], 1.25), div(num_samples, 2))';
Xneg = rand(rng, MvNormal([-1.5, -1.5], 1.25), div(num_samples, 2))';
ypos = ones(div(num_samples, 2));
yneg = -ones(div(num_samples, 2));
nothing #hide
````

````@example svm_primal
# Plot dataset
plot(Xpos[:, 1], Xpos[:, 2], color = :red, st=:scatter, markershape = :rect, label = "positive", xlabel = "x1", ylabel = "x2")
plot!(Xneg[:, 1], Xneg[:, 2], color = :blue, st=:scatter, markershape = :circle, label = "negative")
````

with samples $(x_1, x_2, \ldots, x_m) \in \mathbb{R}^2$ and labels $y_i \in \{-1,1\}$.

## Solving SVM as a QP
We want to compute the weights $w$ and bias term $b$ of the (soft-margin) SVM classifier:

```math
\begin{array}{ll}
    \text{minimize}   & \|w\|^2 + \lambda \sum_{i=1}^m \text{max}(0, 1 - y_i(w^\top x_i  - b)),
\end{array}
```
where $\lambda$ is a hyperparameter. This problem can be solved as a quadratic program.
We can rewrite above problem into an optimisation problem in primal form by introducing the auxiliary slack variables $t_i$:

```math
t_i = \text{max}(0, 1 - y_i(w^T x_i  - b)), \quad t_i \geq 0.
```

This allows us to write the problems in standard QP format:
```math
\begin{array}{ll}
    \text{minimize}   & \|w\|^2 + \lambda \sum_{i=1}^m t_i\\
    \text{subject to} & y_i (w^\top x_i - b) \geq 1 - t_i, \quad \text{for } i = 1,\ldots, m\\
                      & t_i \geq 0, \quad \text{for } i = 1,\ldots, m.
\end{array}
```
Next, we will remove the bias term $b$ by adding an initial feature $x_0 = -1$ to each sample (now: $n = 3$):

````@example svm_primal
X = [-ones(num_samples) [Xpos; Xneg]];
y = [ypos; yneg];
m, n = size(X)
````

## Modelling in JuMP
We can model this problem using `JuMP` and then hand it to `COSMO`:

````@example svm_primal
using JuMP, COSMO
````

````@example svm_primal
λ = 1.0; # hyperparameter
model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "verbose" => true));


@variable(model, w[1:n]);
@variable(model, t[1:m] >= 0.);
@objective(model, Min, w' * w  + λ * ones(m)' * t);
@constraint(model, diagm(0 => y) * X * w .+ t .- 1 .>= 0);
status = JuMP.optimize!(model)
````

The optimal weights $w = [w_0, w_1, w_2]^\top$ (where $w_0 = b$) are:

````@example svm_primal
w_opt = JuMP.value.(w)
````

## Plotting the hyperplane
The separating hyperplane is defined by $w^\top x - b = 0$. To plot the hyperplane, we calculate $x_2$ over a range of $x_1$ values:
```math
x_2 = (-w_1 x_1 - w_0) / w_2, \text{ where } w_0 = b.
```

````@example svm_primal
x1 = -4:0.1:4;
x2 = (-w_opt[2] * x1  .- w_opt[1]) / w_opt[3]
plot!(x1, x2, label = "SVM separator", legend = :topleft)
````

## Modelling with COSMO
The problem can also be solved by transforming it directly into `COSMO`'s problem format.
Define `COSMO``s $x$-variable to be $x=[w, t]^\top$ and choose $P$, $q$, accordingly:

````@example svm_primal
P = blockdiag(spdiagm(0 => ones(n)), spzeros(m, m));
q = [zeros(n); 0.5 * λ * ones(m)];
nothing #hide
````

Next we transform the first constraint $y_i (w^\top x_i - b) \geq 1 - t_i, \quad \text{for } i = 1,\ldots, m$ into
`COSMO`'s constraint format: $Ax + b \in \mathcal{K}$.

````@example svm_primal
A1 = [(spdiagm(0 => y) * X) spdiagm(0 => ones(m))];
b1 = -ones(m);
cs1 = COSMO.Constraint(A1, b1, COSMO.Nonnegatives);
nothing #hide
````

It remains to specify the constraint $t_i \geq 0, \quad \text{for } i = 1,\ldots, m$:

````@example svm_primal
A2 = spdiagm(0 => ones(m));
b2 = zeros(m);
cs2 = COSMO.Constraint(A2, b2, COSMO.Nonnegatives, m+n, n+1:m+n);
nothing #hide
````

Create, assemble and solve the `COSMO.Model`:

````@example svm_primal
model2 = COSMO.Model();
assemble!(model2, P, q, [cs1; cs2]);
result2 = COSMO.optimize!(model2);
w_opt2 = result2.x[1:3];
@test norm(w_opt2 - w_opt, Inf) < 1e-3
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

