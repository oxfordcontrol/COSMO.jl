# # Model Updates
#
# When a pareto-optimal front in portfolio optimisation is computed for different risk-aversion parameters $\gamma$, we are repeatedly solving the same model with only a slight change in the cost function. This allows us to reuse the KKT factorisation from previous solves and recompute the solution very efficiently, see [Algorithm](@ref).
# For more information about solving portfolio optimsiation problems with COSMO see [Portfolio Optimisation](@ref).
# We are planning to solve the following QP for different values of $\gamma$:
# $$
# \begin{array}{ll} \text{minimize} &  x^\top D x + y^\top y - \gamma^{-1} \mu^\top x \\
# \text{subject to} &  y = F^\top x \\
#                   &  1^\top x = d + 1^\top x^0 \\
#                    &  x  \geq 0,
# \end{array}
# $$
# As you can see $\gamma$ enters only in the linear part of the cost function, i.e. our vector `q`. Let's generate some example data:
using LinearAlgebra, SparseArrays, Random, COSMO

rng = Random.MersenneTwister(1)
k = 5; # number of factors
n = k * 10; # number of assets
D = spdiagm(0 => rand(rng, n) .* sqrt(k))
F = sprandn(rng, n, k, 0.5); # factor loading matrix
μ = (3 .+ 9. * rand(rng, n)) / 100. # expected returns between 3% - 12%
d = 1 # we are starting from all cash
x0 = zeros(n);

# We want to solve this problem for the following risk-parameters:
gammas = [ 0.001, 0.01, 0.1,  0.5,  1., 3., 10, 100, 1000]

# Let's solve the problem for the first parameter:
## x = (x, y)
γ = gammas[1]
P = blockdiag(D, spdiagm(0 => ones(k)))
q = [-μ ./ (2 * γ); zeros(k)];

Aeq = [ transpose(F) -diagm(0 => ones(k));
        ones(1, n) zeros(1, k)]
beq = [zeros(k); -1.0]
Aineq = [spdiagm(0 => ones(n)) spzeros(n, k)]

cs1 = COSMO.Constraint(Aeq, beq, COSMO.ZeroSet);
cs2 = COSMO.Constraint(Aineq, zeros(n), COSMO.Nonnegatives);
cosmo_settings = COSMO.Settings(verbose = false)

## solve the problem once
model = COSMO.Model();
assemble!(model, P, q, [cs1; cs2], settings = cosmo_settings);
result = COSMO.optimize!(model);
x_opt = view(result.x, 1:n);
y_opt = view(result.x, n+1:n+k);

## allocate some space for results
risks = zeros(length(gammas));
returns = zeros(length(gammas));
returns[1] = dot(μ, x_opt);
risks[1] = sqrt(dot(y_opt, y_opt));



# After a successful solve, we can use the `update!(model, q = nothing, b = nothing)` function to update our model and re-solve the problem efficiently:


## solve the problem efficiently for the remaining γ using model updates
function solve_repeatedly!(returns, risks, gammas, model, k, n, μ)
        for (i, γ) in enumerate(gammas[2:end])
                q_new = [-μ ./ (2 * γ); zeros(k)];

                ## here we use the update! function to change the model
                update!(model, q = q_new)
                result = COSMO.optimize!(model);

                x_opt = view(result.x, 1:n)
                y_opt = view(result.x, n+1:n+k)
                returns[i + 1] = dot(μ, x_opt)
                risks[i + 1] = sqrt(dot(y_opt, y_opt))
        end
end
solve_repeatedly!(returns, risks, gammas, model, k, n, μ);

# We can now plot the risk-return trade-off curve:
using Plots
Plots.plot(risks, returns, xlabel = "Standard deviation (risk)", ylabel = "Expected return", title = "Risk-return trade-off for efficient portfolios", legend = false)
