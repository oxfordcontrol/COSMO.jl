
# # Support Vector Machine
# We are showing how to solve a support vector machine problem with COSMO (and JuMP).
#
# ## Generating the Dataset
# We want to classify the points in this example dataset with $m = 100$ samples and $n = 2$ features:

using Distributions: MvNormal
using Plots, LinearAlgebra, SparseArrays, Random, Test

#-

## Generate dataset
rng = Random.MersenneTwister(123);
num_samples = 100;
Xpos = rand(rng, MvNormal([1.5, 1.5], 1.25), div(num_samples, 2))';
Xneg = rand(rng, MvNormal([-1.5, -1.5], 1.25), div(num_samples, 2))';
ypos = ones(div(num_samples, 2));
yneg = -ones(div(num_samples, 2));

#-
## Plot dataset
plot(Xpos[:, 1], Xpos[:, 2], color = :red, st=:scatter, markershape = :rect, label = "positive", xlabel = "x1", ylabel = "x2")
plot!(Xneg[:, 1], Xneg[:, 2], color = :blue, st=:scatter, markershape = :circle, label = "negative")

# with samples $(x_1, x_2, \ldots, x_m) \in \mathbb{R}^2$ and labels $y_i \in \{-1,1\}$.
#
# ## Solving SVM as a QP
# We want to compute the weights $w$ and bias term $b$ of the (soft-margin) SVM classifier:
#
# $$
# \begin{array}{ll}
#     \text{minimize}   & \|w\|^2 + \lambda \sum_{i=1}^m \text{max}(0, 1 - y_i(w^\top x_i  - b)),
# \end{array}
# $$
# where $\lambda$ is a hyperparameter. This problem can be solved as a quadratic program.
# We can rewrite above problem into an optimisation problem in primal form by introducing the auxiliary slack variables $t_i$:
#
# $$
# t_i = \text{max}(0, 1 - y_i(w^T x_i  - b)), \quad t_i \geq 0.
# $$

# This allows us to write the problems in standard QP format:
# $$
# \begin{array}{ll}
#     \text{minimize}   & \|w\|^2 + \lambda \sum_{i=1}^m t_i\\
#     \text{subject to} & y_i (w^\top x_i - b) \geq 1 - t_i, \quad \text{for } i = 1,\ldots, m\\
#                       & t_i \geq 0, \quad \text{for } i = 1,\ldots, m.
# \end{array}
# $$
# Next, we will remove the bias term $b$ by adding an initial feature $x_0 = -1$ to each sample (now: $n = 3$):

X = [-ones(num_samples) [Xpos; Xneg]];
y = [ypos; yneg];
m, n = size(X)

# ## Modelling in JuMP
# We can model this problem using `JuMP` and then hand it to `COSMO`:
using JuMP, COSMO
#-
λ = 1.0; # hyperparameter
model = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose=true));


@variable(model, w[1:n]);
@variable(model, t[1:m] >= 0.);
@objective(model, Min, w' * w  + λ * ones(m)' * t);
@constraint(model, diagm(0 => y) * X * w .+ t .- 1 .>= 0);
status = JuMP.optimize!(model)
# The optimal weights $w = [w_0, w_1, w_2]^\top$ (where $w_0 = b$) are:
w_opt = JuMP.value.(w)


# ## Plotting the hyperplane
# The separating hyperplane is defined by $w^\top x - b = 0$. To plot the hyperplane, we calculate $x_2$ over a range of $x_1$ values:
# $$
# x_2 = (-w_1 x_1 - w_0) / w_2, \text{ where } w_0 = b.
# $$
x1 = -4:0.1:4;
x2 = (-w_opt[2] * x1  .- w_opt[1]) / w_opt[3]
plot!(x1, x2, label = "SVM separator", legend = :topleft)

# ## Modelling with COSMO
# The problem can also be solved by transforming it directly into `COSMO`'s problem format.
# Define `COSMO``s $x$-variable to be $x=[w, t]^\top$ and choose $P$, $q$, accordingly:
P = blockdiag(spdiagm(0 => ones(n)), spzeros(m, m));
q = [zeros(n); 0.5 * λ * ones(m)];

# Next we transform the first constraint $y_i (w^\top x_i - b) \geq 1 - t_i, \quad \text{for } i = 1,\ldots, m$ into
# `COSMO`'s constraint format: $Ax + b \in \mathcal{K}$.
A1 = [(spdiagm(0 => y) * X) spdiagm(0 => ones(m))];
b1 = -ones(m);
cs1 = COSMO.Constraint(A1, b1, COSMO.Nonnegatives);

# It remains to specify the constraint $t_i \geq 0, \quad \text{for } i = 1,\ldots, m$:
A2 = spdiagm(0 => ones(m));
b2 = zeros(m);
cs2 = COSMO.Constraint(A2, b2, COSMO.Nonnegatives, m+n, n+1:m+n);

# Create, assemble and solve the `COSMO.Model`:
model2 = COSMO.Model();
assemble!(model2, P, q, [cs1; cs2]);
result2 = COSMO.optimize!(model2);
w_opt2 = result2.x[1:3];
@test norm(w_opt2 - w_opt, Inf) < 1e-3
