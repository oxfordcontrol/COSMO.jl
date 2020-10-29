# # Portfolio Optimisation
#
# We consider a single-period Markowitz portfolio optimisation example.
# Assume that we have a portfolio with $n$ assets at the beginning of time period $t$. Given some forecasts on risks and expected returns we try to find the optimal trade vector that rebalances the portfolio to achieve a good balance between expected risk (variance) $x^\top \Sigma x$ and returns $\mu^\top x$. In it's most simple form we want to solve:
# $$
# \begin{array}{ll} \text{maximize} &  \mu^\top x - \gamma (x^\top \Sigma x)\\
# \text{subject to} &  1^\top x = d + 1^\top x^0 \\
#                    &  x  \geq 0,
# \end{array}
# $$
# with variable $x \in \mathbf{R}^n$, $\mu$ forecasted (expected) returns, $\gamma > 0 $ risk aversion parameter. $x^0_i$ represents the initial investment in asset $i$ and $d$ represents the cash reserve. Consequently, the equality constraint tells us that the sum of the new allocation vector $x$ has to equal the initial allocation plus the cash reserve. Furthermore, the covariance matrix of our risk model is given by $\Sigma \in \mathbf{S}_+^n$. Here we assume a factor risk model, i.e. we can write $\Sigma$ as $\Sigma = D + F F^\top$ where $D$ is diagonal and the factor matrix $F$ has a lower rank $k < n$. This approach allows us to reduce the number of nonzeros in the problem. Furthermore, note that we don't consider shortselling in this example. Let's generate some problem data:
#-
using LinearAlgebra, SparseArrays, Random, COSMO, JuMP, Test

## generate the data
rng = Random.MersenneTwister(1)
k = 5; # number of factors
n = k * 10; # number of assets
D = spdiagm(0 => rand(rng, n) .* sqrt(k))
F = sprandn(rng, n, k, 0.5); # factor loading matrix
μ = (3 .+ 9. * rand(rng, n)) / 100. # expected returns between 3% - 12%
γ = 1.0; # risk aversion parameter
d = 1 # we are starting from all cash
x0 = zeros(n);

# We can now write the problem as a QP:
# $$
# \begin{array}{ll} \text{minimize} &  x^\top D x + y^\top y - \gamma^{-1} \mu^\top x \\
# \text{subject to} &  y = F^\top x \\
#                   &  1^\top x = d + 1^\top x^0 \\
#                    &  x  \geq 0.
# \end{array}
# $$
# Before considering other effects, let's create the model in JuMP and solve it using COSMO:
#-

model = JuMP.Model(with_optimizer(COSMO.Optimizer));
@variable(model, x[1:n]);
@variable(model, y[1:k]);
@objective(model, Min,  x' * D * x + y' * y - 1/γ * μ' * x);
@constraint(model, y .== F' * x);
@constraint(model, sum(x) == d + sum(x0));
@constraint(model, x .>= 0);
JuMP.optimize!(model)

# After solving the problem, we can calculate the expected return and risk $\sigma= \sqrt{x^{* \top} \Sigma x^*}$:
x_opt = JuMP.value.(x);
y_opt = JuMP.value.(y);
expected_return_basic = dot(μ, x_opt)

#-
expected_risk_basic = sqrt(dot(y_opt, y_opt))

# ## Using standard deviation in the model
# It is pointed out in \[1\] that above problem formulation can lead to numerical problems, e.g. if $\Sigma$ is not strictly positive semidefinite. Another option is to formulate the risk constraint in terms of the standard deviation $\|M^\top x \|$ where $M M^\top = D + F F^\top$ and bound it using a second-order cone constraint:
# $$
# \begin{array}{ll} \text{minimize} &  - \mu^\top x \\
# \text{subject to} &    \|M^\top x\| \leq \gamma \\
#                   &  1^\top x = d + 1^\top x^0 \\
#                    &  x  \geq 0.
# \end{array}
# $$

Mt = [D.^0.5; F']
model = JuMP.Model(with_optimizer(COSMO.Optimizer));
@variable(model, x[1:n]);
@objective(model, Min, - μ' * x);
@constraint(model,  [γ; Mt * x] in SecondOrderCone()); # ||M'x|| <= γ
@constraint(model, sum(x) == d + sum(x0));
@constraint(model, x .>= 0);
JuMP.optimize!(model)

# Note that the result is different from the example above because $\gamma$ scales the problem in a different way. Here it can be seen as an upper bound on the standard deviation of the portfolio.
x_opt = JuMP.value.(x);
expected_return = dot(μ, x_opt)
# Let us verify that the bound holds:
@test norm(Mt * x_opt) <= γ


# ## Pareto-optimal front
# The above portfolio optimisation approach yields the optimal expected return for a given level of risk. The result is obviously impacted by the risk aversion $\gamma$ parameter. To visualise the trade-off and present the investor with an efficient Pareto optimal portfolio for their risk appetite we can compute the optimal portfolio for many choices of $\gamma$ and plot the corresponding risk-return trade-off curve.

gammas = [ 0.001, 0.01, 0.1,  0.5,  1., 3., 10, 100, 1000]
risks = zeros(length(gammas))
returns = zeros(length(gammas))
model = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose = false));
@variable(model, x[1:n]);
@variable(model, y[1:k]);
@objective(model, Min,  x' * D * x + y' * y - 1/γ * μ' * x);
@constraint(model, y .== F' * x);
@constraint(model, sum(x) == d + sum(x0));
@constraint(model, x .>= 0);

## solve the same problem for different values of γ
for (k, gamma) in enumerate(gammas)
    coeff = - 1/gamma * μ
    JuMP.set_objective_coefficient.(model, x, coeff)
    JuMP.optimize!(model)
    local x_opt = JuMP.value.(x);
    local y_opt = JuMP.value.(y);
    returns[k] = dot(μ, x_opt)
    risks[k] = sqrt(dot(y_opt, y_opt))
end
# We can now plot the risk-return trade-off curve:
using Plots
Plots.plot(risks, returns, xlabel = "Standard deviation (risk)", ylabel = "Expected return", title = "Risk-return trade-off for efficient portfolios", legend = false)

# !!! note
#     When the model is updated in `JuMP` as above the `JuMP.model` is copied in full to `COSMO`. We are trying to improve the interface with respect to model updates in the future. Until then you can use [Model Updates](@ref) in `COSMO`s native interface.



# ## Transaction costs
# In the model above we assume that trading the assets is free and does not impact the market. However, this is clearly not the case in reality. To make the example more realistic consider the following cost $c_j$ associated with the trade $δ_j = x_j - x_j^0$:
# $$
# c_j(\delta_j) = a_j |\delta_j| + b_j |\delta_j|^{3/2},
# $$
# where the first term models the bid-ask spread and broker fees for asset $j$. The second term models the impact on the market that our trade has. This is obviously only a factor if the volume of our trade is significant. The constant $b_j$ is a function of the total volume traded in the considered time periode and the price volatility of the asset and has to be estimated by the trader. To make this example simple we consider the same coefficients $a$ and $b$ for every asset. The $|\delta_j|^{3/2}$ term can be easily modeled using a power cone constraint $\mathcal{K}_{pow} = \{(x, y, z) \mid x^\alpha y^{(1-\alpha)} \geq |z|, x \geq 0, y \geq 0, 0 \leq \alpha \leq 1 \}$. In fact this can be used to model any market impact function with exponent greater than 1.
# We can write the total transaction cost $a^\top s + b^\top t$ where $s_j$ bounds the absolute value of $\delta_j$ and $t_{j}$ is used to bound the term $|x_j - x_j^0|^{3/2} \leq t_{j}$ using a power cone formulation: $(t_{j}, 1, x_j - x_j^0) \in \mathcal{K}_{pow}(2/3)$.

#-
a = 1e-3
b = 1e-1
γ = 1.0;
model = JuMP.Model(with_optimizer(COSMO.Optimizer, eps_abs = 1e-5, eps_rel = 1e-5));
@variable(model, x[1:n]);
@variable(model, y[1:k]);
@variable(model, t[1:n]);
@variable(model, s[1:n]);
@objective(model, Min, x' * D * x + y' * y - 1/γ * μ' * x);
@constraint(model, y .== F' * x);
@constraint(model, x .>= 0);

## transaction costs
@constraint(model, sum(x) + a * sum(s) + b * sum(t) == d + sum(x0) );
@constraint(model, [i = 1:n], x[i] - x0[i] <= s[i]); # model the absolute value with slack variable s
@constraint(model, [i = 1:n], x0[i] - x[i] <= s[i]);
@constraint(model, [i = 1:n], [t[i], 1, x[i] - x0[i]] in MOI.PowerCone(2/3));
JuMP.optimize!(model)
# Let's look at the expected return and the total transaction cost:

x_opt = JuMP.value.(x);
y_opt = JuMP.value.(y);
s_opt = JuMP.value.(s);
t_opt = JuMP.value.(t);
expected_return = dot(μ, x_opt)
#-
expected_risk = dot(y_opt, y_opt)
#-
transaction_cost = a * sum(s_opt) + b * sum( t_opt)


# ## References
# [1] [Mosek Case Studies](https://docs.mosek.com/9.2/pythonfusion/case-studies-portfolio.html)
