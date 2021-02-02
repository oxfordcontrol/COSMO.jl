# # Lovász Theta Function
#
# The Lovász theta function is an important concept in graph theory. Consider an undirected graph $G = (V,E)$ with vertex set $V = \{1,\dots,n\}$ and edge set $E$ with  $E_{ij} = 1$ if there exists an edge between the vertices $i$ and $j$. A stable set (or independent set) is a subset of $V$ such that the induced subgraph does not contain edges. The *stability number* $\alpha(G)$ of the graph is equal to the cardinality of the largest stable set. It is closely related to the *Shannon capacity* $\Theta(G)$ that models the amount of information that a noisy communication channel can carry if certain signal values can be confused with each other. Unfortunately, the determination of $\alpha(G)$ is an NP-hard problem and the computational complexity to determine $\Theta(G)$ remains unknown. The Lovász theta function $\vartheta(G)$ was introduced in \[1\], can be computed in polynomial time, and represents an upper bound on both the stability number and the Shannon capacity:
# $$
# \alpha(G) \leq \Theta(G) \leq \vartheta(G).
# $$
# The value of $\vartheta(G)$ can be determined by computing the optimal value $p^*$ of the following SDP:
# $$
# \begin{array}{ll} \text{maximize} &   \text{Tr}(JX)\\
# \text{subject to} &  \text{Tr}(X) = 1, \\
#                   &  X_{ij} = 0, \quad (i, j) \in E, \\
#                   &  X \succeq 0,
# \end{array}
# $$
# with matrix variable $X$, matrix $J \in \mathbf{R}^{n \times n}$ of all ones and $\text{Tr}()$ denoting the matrix trace.
# In this simple example we will compute the value of the Lovász theta function for the *Petersen Graph* (which is known to be $4$).

#-
using LinearAlgebra, COSMO, JuMP, Plots, GraphRecipes

## let's define the Petersen graph
n = 10
E = zeros(n, n)
E[1, 2] = E[1, 5] = E[1, 6] = 1.
E[2, 3] = E[2, 7]  = 1.
E[3, 4] = E[3, 8]  = 1.
E[4, 5] = E[4, 9] = 1.
E[5, 10] = 1.
E[6, 8] = E[6, 9] = 1.
E[7, 9] = E[7, 10] = 1.
E[8, 10] = 1.

## plot the graph
ri = 1.
ro = 2.
coordinates = []
for θ = 90:-72:-198
  push!(coordinates, [ro * cosd(θ), ro * sind(θ)])
end
for θ = 90:-72:-198
  push!(coordinates, [ri * cosd(θ), ri * sind(θ)])
end
graphplot(E, names = 1:n,  x = getindex.(coordinates, 1), y = getindex.(coordinates, 2), fontsize = 10, nodesize = 1, nodeshape =:circle, curvature = 0.)

# Let's solve the SDP with COSMO and JuMP:
#-
model = JuMP.Model(COSMO.Optimizer);

@variable(model, X[1:n, 1:n], PSD)
x = vec(X)
@objective(model, Max, sum(x))
@constraint(model, tr(X)== 1.)
for j = 1:n
  for i = 1:j-1
    if E[i, j] == 1.
      @constraint(model, X[i, j] == 0)
    end
  end
end
status = JuMP.optimize!(model)

# The optimal objective is given by:
#-
JuMP.objective_value(model)
# Which is the correct known value for the Petersen Graph.
#
# ## References
# [1] Lovász - On the Shannon Capacity of a Graph, IEEE Transactions on Information Theory (1979)
