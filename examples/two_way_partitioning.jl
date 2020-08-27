# # Relaxed Two-Way Partitioning Problem
#
# We consider the (nonconvex) two-way partitioning problem from Boyd and Vandenberghe (2004), p.219 \[1\]:
# $$
# \begin{array}{ll} \text{minimize} &  x^\top W x \\
# \text{subject to} &  x_i^2 = 1, \quad i=1,\dots, n,
# \end{array}
# $$
# with $x \in \mathbf{R}^n$ and $W \in \mathbf{S}^n$. The problem can be interpreted as finding a partition of the points $i = 1,\dots, n$ into two sets, where the cost of two points in the same set is $W_{ij}$ and $-W_{ij}$ otherwise. Brute forcing the solution $x^*$ takes $2^n$ attempts and becomes quickly intractable, e.g. for $n \geq 30$. However, a lower-bound for this problem can be computed by solving the convex dual of the problem:
# $$
# \begin{array}{ll} \text{maximize} &  -\sum_i \nu \\
# \text{subject to} &  W + \text{diag}(\nu) \in \mathbf{S}_+^n.
# \end{array}
# $$
# Solving this problem with optimal objective $d^*$ yields a lower bound as least as good as a lower bound based on the minimum eigenvalue $d^* \geq \lambda_{\text{min}}(W)$.
# Let's set up the problem for $n = 20$ and with a specially structured (banded) $W$:

#-
using LinearAlgebra, Random, COSMO, IterTools, JuMP

rng = Random.MersenneTwister(212313);
n = 20;

W = diagm(0 => randn(rng, n));
W += diagm(1 => randn(rng, n-1));
W = Symmetric(W)

# As you can see, the matrix $W$ imposes a structure on the LMI-constraint. Accordingly, COSMO will be able to use chordal decomposition to decompose the LMI constraint into multiple smaller constraints. This can make a significant difference in the performance of the algorithm for large $n$.

#-

model = JuMP.Model(with_optimizer(COSMO.Optimizer));
@variable(model, ν[1:n]);
@objective(model, Max, -ones(n)' * ν )
@constraint(model, Symmetric(W + diagm(ν) )  in JuMP.PSDCone());
JuMP.optimize!(model)

# Looking at the solver output you can see how the PSD constraint was decomposed into 19 PSD constraints. Let's look at the lower bound:
JuMP.objective_value(model)

# As $n$ is small, we can verify our result by finding the optimal solution by trying out all possible combinations:

function brute_force_optimisation(W, n)
   opt_obj = Inf
   opt_x = Inf * ones(n)

   for xt in Iterators.product([[1.; -1.] for k = 1:n]...)
      x = [i for i in xt]
      obj_val = x' * W * x
      if obj_val < opt_obj
         opt_obj = obj_val
         opt_x = x
      end
   end
   return opt_obj, opt_x
end
opt_obj, opt_x = brute_force_optimisation(W, n)
# ## References
# [1] Boyd and Vandenberghe - Convex Optimization, Cambridge University Press (2004)
