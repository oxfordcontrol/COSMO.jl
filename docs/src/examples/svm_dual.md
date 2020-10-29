The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
```@meta
EditURL = "<unknown>/../examples/svm_dual.jl"
```

# Support Vector Machine (dual)
We are showing how to solve the dual support vector machine formulation with the kernel trick. We solve the resulting QP with COSMO (and JuMP).

## Loading the dataset
Let's take a look at the dataset for social network ads. We want to build a classifier based on the age and salary of a customer and then decide whether they will purchase based on the ad shown.

```@example svm_dual
using CSV, DataFrames, Plots, LinearAlgebra, StatsBase
```

```@example svm_dual
# load data set into a dataframe
df = CSV.read("./data/social_network_ads.csv")
Xd = convert(Array{Float64, 2}, [df.Age df.EstimatedSalary])
y = df.Purchased[:]
y[y .== 0] .= -1 # -1: not purchased, 1: purchased
```

```@example svm_dual
# Plot dataset
scatter(Xd[y .== -1, 1], Xd[y .== -1, 2],  color = :red, xlabel = "Age", ylabel = "Esimated Salary", label = "not purchased", title = "Social Network Ad data", legend=:topleft)
scatter!(Xd[y .== 1, 1], Xd[y .== 1, 2],  color = :green, label = "purchased")
```

with samples $(x_1, x_2, \ldots, x_m) \in \mathbb{R}^2$ and labels $y_i \in \{-1,1\}$. We can see that the data will not be linearly separable. We are therefore interested in using a kernel function to transform the feature space.
Before we continue we will standardize the dataset to have 0 mean and unit variance.

```@example svm_dual
dt = fit(ZScoreTransform, Xd, dims=1)
Xstd = StatsBase.transform(dt, Xd)
n, m = size(Xstd)
```

## Deriving the dual SVM formulation
[put all the derivations here]

define kernel function

```@example svm_dual
RBF(x1::AbstractVector{T}, x2::AbstractVector{T}, σ::T) where {T <: AbstractFloat} = exp(-norm(x1 - x2)^2 / (2 * σ^2))
```

## Modelling in JuMP
We can model this problem using `JuMP` and then hand it to `COSMO`:

```@example svm_dual
using JuMP, COSMO
```

```@example svm_dual
C = 10.; # constraint-violation penalty term (in primal problem)
σ = 1.
```

compute kernel matrix

```@example svm_dual
K = zeros(n, n)
for c = 1:n, r = 1:c
    x1 = @view Xstd[c, :]
    x2 = @view Xstd[r, :]
    K[r, c] = RBF(x1, x2, σ)
end
K = Symmetric(K, :U)

model = JuMP.Model(with_optimizer(COSMO.Optimizer, verbose=true));

@variable(model, α[1:n]);
@objective(model, Max, sum(α) - 0.5 * (α .* y)' * K * (α .* y));
@constraint(model, sum(α .* y) == 0);
@constraint(model, 0. .<= α .<= C);
status = JuMP.optimize!(model)
```

The optimal alphas are:

```@example svm_dual
α_opt = JuMP.value.(α)
```

find the alphas that are nonzero corresponding to the support vectors

```@example svm_dual
SV = α_opt .>= 1e-3
SV = findall(x -> x == 1, SV)

α_b = (1e-3 .<= α_opt .<= C)
α_b = findall(x -> x == 1, α_b)

function compute_bias(Xstd, y, SV, σ, α_opt, α_b)
    b = 0.
    for ind in α_b
        xi =  @view Xstd[ind, :]
        wTx = 0.
        for SV_ind in SV
            xj = @view Xstd[SV_ind, :]
            wTx += α_opt[SV_ind] * y[SV_ind] * RBF(xi, xj, σ)
        end
        b += (y[ind] - wTx)
    end
    b /= length(α_b)
    return b
end
b = compute_bias(Xstd, y, SV, σ, α_opt, α_b)
```

## Plotting the hyperplane
The separating hyperplane is defined by $w^\top x - b = 0$. To plot the hyperplane, we calculate $x_2$ over a range of $x_1$ values:
```math
x_2 = (-w_1 x_1 - w_0) / w_2, \text{ where } w_0 = b.
```

```@example svm_dual
function eval_svm(u::Vector{T}, α::Vector{T}, SV::Vector{Int64}, b::T, Xstd::Matrix{T}, y::Vector{Int64}, σ::T) where {T <: AbstractFloat}
    z = b
    for ind in SV
        xi = @view Xstd[ind, :]
        z += α[ind] * y[ind] * RBF(xi, u, σ)
    end
    return z
end

u = -3:0.01:3;
v = -3:0.01:3;
z = zeros(length(u), length(v));
for i = 1:length(u), j = 1:length(v)
    val = eval_svm([u[i]; v[i]], α_opt, SV, b, Xstd, y, σ)
    z[i, j] = sign(val);#sign(val)
end
U = [i for i in u, j in v]
V = [j for i in u, j in v]
```

plot the decision boundary

```@example svm_dual
PyPlot.contour(u, v, z', levels = [1.])
```

contour(u, v, sign(eval_svm(u, vz', f= true)# c = :black, linewidth = 2)
scatter!(Xstd[y .== -1, 1], Xstd[y .== -1, 2],  color = :red, xlabel = "Age", ylabel = "Esimated Salary", label = "not purchased", title = "Social Network Ad data", legend=:topleft)
scatter!(Xstd[y .== 1, 1], Xstd[y .== 1, 2],  color = :green, label = "purchased")

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

