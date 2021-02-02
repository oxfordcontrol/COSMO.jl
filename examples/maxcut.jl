# # Maximum Cut Problem
#
# We are interested in solving an approximation to the maximum cut problem using semidefinite programming.
# Consider the graph $G(V,E)$ with weights $w_{ij}$ shown below:

## create the weighted adjacency matrix of the graph and plot the graph
using LinearAlgebra, Plots, GraphRecipes
n = 4;
W = zeros(n, n);
W[1, 2] = 1; W[1, 4] = 8;
W[2, 3] = 2; W[2, 4] = 10;
W[3, 4] = 6;
W = Symmetric(W)
graphplot(W, names = 1:n, edgelabel = W,  x = [0; 1; 1; 0], y = [1; 1; 0; 0], fontsize = 12, nodeshape =:circle)

# The maximum cut problem tries to find a cut or a partition of the graph's vertices into two complementary sets $S$ and $\bar{S}$ such that the total weight of the edges between the two sets is maximized. For this small graph the problem is trivial. The optimal solution is $S = \{1,2,3 \}$ and $\bar{S}=\{4\}$.
# Formally, this problem can be written as a mixed-integer optimisation problem:
# $$
# \begin{array}{ll} \text{maximize} &  \frac{1}{2} \sum_{i < j} w_{ij}(1 - y_i y_j)\\
# \text{subject to} &  y_i \in \, \{-1, 1 \}, \quad \forall i \in V,
# \end{array}
# $$
# where $y_i = 1$ indicates that $v_i \in S$ and $y_i = -1$ indicates that $v_i \in \bar{S}$. This problem is of interest in the field of integrated circuit design, where one tries to minimize the number of cross-layer connections in a circuit under layout constraints.

# For more complicated graphs this problem quickly becomes hard to solve to optimality and in fact is known to be NP-hard. For this example we are interested in the randomized approximation algorithm that relaxes the problem to an SDP and was devised in Goemans and Williamson (1995) \[1\] (see for more details).

# The approach can be divided into three steps:
# 1. Solve a relaxed SDP to obtain $Y^*$
# 2. Recover an approximate solution $V$ via a Cholesky factorisation $Y^* = V^\top V$
# 3. Round the approximate solution using a random vector $r$ from the unit sphere.
# The authors showed that this approach guarantees a solution of at least $0.87856$ times the optimal solution.

# ## Solving the primal SDP
# Before we formulate the SDP, let's compute the Laplacian matrix $L$:
using LinearAlgebra, Random, SparseArrays, COSMO, JuMP

rng = Random.MersenneTwister(1)

## Compute the Laplacian matrix of the graph
L = diagm(0 => W * ones(n)) - W;

# Given the Laplacian matrix $L$, we are looking for the optimal $Y^*$ that solves the following SDP:
# $$
# \begin{array}{ll} \text{maximize} &  \frac{1}{4} \langle L, Y \rangle \\
# \text{subject to} & Y_{ii} = 1 \quad \text{for } i = 1,\dots,n\\
#                    & Y \in \mathbf{S}_{+}^n.
# \end{array}
# $$
# Notice that the solution $Y^*$ can be viewed as a correlation matrix. Let's solve the problem using COSMO and JuMP.

#-
model = JuMP.Model(COSMO.Optimizer);
@variable(model, Y[1:n, 1:n], PSD);
@objective(model, Max, 1 / 4 * dot(L, Y));
@constraint(model, [i = 1:n], Y[i, i] == 1.);
JuMP.optimize!(model)

#-
Yopt = JuMP.value.(Y);
obj_val = JuMP.objective_value(model)

# ## Solving the dual SDP
# Notice that the decision matrix $Y$ is generally dense (as correctly classified in the solver output above). Therefore, we won't be able to utilize COSMO's chordal decomposition features. However, assuming strong duality, it turns out that we can also solve the dual problem, which is given by:
# $$
# \begin{array}{ll} \text{minimize} &  \sum_i \gamma_i \\
# \text{subject to} & \text{diag}(\gamma) - \frac{1}{4} L = S \\
#                   & S \in \mathbf{S}_{+}^n.
# \end{array}
# $$
# As you can see, the matrix $S$ is constrained to have zeros in places specified by the graph (Laplacian). Therefore, COSMO can try to decompose this SDP (see solver output) and speed up its algorithm:

#-
model_dual = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "complete_dual" => true));
@variable(model_dual, γ[1:n]);
@objective(model_dual, Min,  sum(γ));
@constraint(model_dual, lmi, Symmetric(-1/4 * L + diagm(γ)) in JuMP.PSDCone());
JuMP.optimize!(model_dual)
obj_val = JuMP.objective_value(model_dual)

# The primal variable $Y^*$ can be recovered from the dual variable associated with the LMI-constraint:
#-
Yopt = dual(lmi)
# To get the correct positive semidefinite dual variable, we have to enable PSD completion of the dual variable in COSMO. Now, that we have a solution we can perform the remaining steps in the approximation algorithm.

# ## Cholesky factorisation of Y
# Compute the Cholesky factorisation of $Y = V^\top V$ to find the unit vectors $v_1, \dots, v_n$.

#-
factor = cholesky(Yopt, Val(true); check = false);
V = Matrix((factor.P * factor.L)');
## normalize columns
for i in 1:n
    V[:, i] ./= norm(V[:, i]);
end

# ## Rounding the approximate solution
# It remains to round the unit vectors $v_i$ using a random vector $r$ with each component drawn from $\mathcal{U}(0, 1)$, to obtain the $y_i$'s:

#-
r = rand(rng, n)
r /= norm(r, 2)
y = ones(n)
for i in 1:n
  dot(r, V[:, i]) <= 0 && (y[i] = -1)
end
y
# For larger graphs this rounding step could be repeated several times to improve the rounding.
# As expected $S = \{1, 2, 3\}$ and $\bar{S}=\{ 4\}$.
# ## References
# [1] Goemans and Williamson - Improved Approximation Algorithms for Maximum Cut and Satisfiability Problems Using Semidefinite Programs (1995)

# [2] Post-processing code steps 2 and 3 from this [JuMP example](https://github.com/jump-dev/JuMP.jl/blob/master/examples/max_cut_sdp.jl)
