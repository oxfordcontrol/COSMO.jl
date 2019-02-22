# Method
This section describes COSMO's underlying ADMM algorithm and how the user can use the settings to adjust this algorithm. For a more detailed explanation take a look at the associated publication in [Citing COSMO](@ref).

COSMO solves problems with quadratic objective function and a number of conic constraints in the following form:
```math
\begin{array}{ll} \mbox{minimize} & \textstyle{\frac{1}{2}}x^\top Px + q^\top x\\ \mbox{subject to} & Ax + s = b \\ & s \in \mathcal{K}, \end{array}
```

with primal decision variable ``x \in \mathbb{R}^n``, primal slack variable ``s \in \mathbb{R}^m``. The objective function is defined by positive semidefinite matrix ``P=P^\top \succeq 0`` and vector ``q \in \mathbb{R}^n``. The constraints are defined by matrix ``A \in \mathbb{R}^{m \times n}``, vector ``b \in \mathbb{R}^m`` and a non-empty, closed convex set ``\mathcal{K}``. The convex set itself can be a Cartesian product of convex sets in the form:
```math
  \mathcal{K} = \mathcal{K}_1^{m_1} \times \mathcal{K}_2^{m_2} \times \cdots \times \mathcal{K}_N^{m_N}.
```
Accordingly, by an appropriate choice of convex sets one can represent any LP, QP, SOCP or SDP.

## Dual problem

The dual problem of the optimisation problem above is given by:
```math
\begin{alignat*}{2}
&\text{maximize}   & \quad   & -\textstyle{\frac{1}{2}}x^\top Px - b^\top  y - \text{sup}_{s \in \mathcal{K}}(-y^\top s ) \\
&\text{subject to} & & Px + A^\top  y = -q\\
&            & & y \in (\mathcal{K}^\infty)^* ,
\end{alignat*}
```
with dual variable ``y \in \mathbb{R}^m``.


## Algorithm
The algorithm considers a slightly transformed problem. By introducing two dummy variables ``\tilde{x} = x`` and ``\tilde{s} = s`` we can rewrite the original problem:
```math
\begin{alignat*}{2}
&\text{minimize}   &  ~   & \textstyle{\frac{1}{2}}\tilde{x}^\top P \tilde{x} + q^\top \tilde{x} + \mathcal{I}_{Ax+s=b}(\tilde{x},\tilde{s}) + \mathcal{I}_{\mathcal{K}}(s)\\
&\text{subject to} & & (\tilde{x},\tilde{s}) = (x,s),
\end{alignat*}
```
 where indicator functions ``\mathcal{I}`` were used to move the constraints into the objective function. The resulting problem is now in the right format to apply the alternating direction method of multipliers (ADMM). To apply ADMM we first find the augmented Lagrangian ``L``:
```math
  L(x,s,\tilde{x},\tilde{s},\lambda,y) = \textstyle{\frac{1}{2}}\tilde{x}^\top P\tilde{x} + q^\top \tilde{x} + \mathcal{I}_{Ax+s=b}(\tilde{x},\tilde{s}) + \mathcal{I}_{\mathcal{K}}(s) + \frac{\sigma}{2} \|\tilde{x} - x + \textstyle{\frac{1}{\sigma}} λ \|_2^2 + \frac{\rho}{2} \| \tilde{s} - s + \textstyle{\frac{1}{\rho}} y \|_2^2.
```
 Minimizing the Lagrangian in an alternating fashion with respect to the two variable pairs ``(\tilde{x},\tilde{s})`` and ``(x,s)`` yields the following algorithm steps:
```math
\begin{align*}
    ( \tilde{x}^{k+1},\tilde{s}^{k+1})  &\rightarrow \underset{\tilde{x},\tilde{s}}{\text{argmin }}  L\left( \tilde{x},\tilde{s},x^k,s^k,y^k \right)\\
    x^{k+1} &\leftarrow \tilde{x}^{k+1}  \label{eq:ADMM1}\\
    s^{k+1} &\leftarrow \underset{s}{\text{argmin }}\frac{\rho}{2} \|  \tilde{s}^{k+1}  - s+\textstyle{\frac{1}{\rho}}y^k \|_2^2 + I_{\mathcal{K}}(s), \\
    y^{k+1} &\leftarrow y^k + \rho \left( \tilde{s}^{k+1} -s^{k+1} \right).
\end{align*}
```
By the construction of the ADMM method those iterates are converging to the global solution. These steps are executed in a loop until convergence. Two important parameters are the ADMM steps sizes ``\rho`` (`Settings.rho`) and ``\sigma`` (`Settings.sigma`) which can be adjusted via the solver settings.

The two most important steps of the algorithm happen in the first and third line. The evaluation of the first line turns out to be an equality constrained quadratic program. We get a solution for ``( \tilde{x}^{k+1},\tilde{s}^{k+1})`` at every iteration by solving the following linear system:
```math
\begin{align*}
\begin{bmatrix}
P + \sigma I & A^\top \\A &- \frac{1}{\rho}I
    \end{bmatrix}\begin{bmatrix}\tilde{x}^{k+1} \\ \nu^{k+1}\end{bmatrix}&= \begin{bmatrix}-q+\sigma x^k \\b-s^k+\frac{1}{\rho}y^k\end{bmatrix}\\
\tilde{s}^{k+1} &= s^k - \frac{1}{\rho}\left(\nu^{k+1} + y^{k}\right).
\end{align*}
```
Fortunately, the left hand matrix doesn't change, which is why COSMO only has to factor the matrix once at the start of the algorithm.

The second important step in the algorithm is the update equation for ``s^{k+1}`` which can be interpreted as a projection onto the constraint set ``\mathcal{K}``:
```math
s^{k+1} = \Pi_{\mathcal{K}}\left( \tilde{s}^{k+1} + \frac{1}{\rho}y^k\right).
```
The computational cost of this projection is highly dependent on the constraints of the problem. While projections onto the zero set or the nonnegative orthant are inexpensive, projections onto the positive semidefinite cone of order ``N`` involve an eigen-decomposition. Since methods for eigen-decompositions have a complexity of ``\mathcal{O}(N^3)`` the projection can become the computationally most expensive operation of the algorithm.

## Scaling
The convergence of ADMM-based algorithms depends on the relative scaling of the problem data. Especially to improve the convergence of badly scaled problems, COSMO tries to rescale the data in a preprocessing step.

We rescale the equality constraints with diagonal positive semidefinite matrices ``D`` and ``E``. The scaled problem is given by:
```math
\begin{alignat*}{2}
\label{eq:scaled}
&\text{minimize}   & \quad   & \textstyle{\frac{1}{2}} \hat{x}^\top \hat{P} \hat{x} + \hat{q}^\top \hat{x}\\
&\text{subject to} & & \hat{A} \hat{x} + \hat{s}  = \hat{b},  \\
&            & & \hat{s} \in E\mathcal{K},
\end{alignat*}
```
with scaled problem data
```math
  \hat{P}=DPD, \quad \hat{q}=Dq,  \quad\hat{A}=EAD, \quad \hat{b}=Eb,
```
and the scaled convex cone ``E\mathcal{K} = \{Ev \in \mathbb{R}^m \mid v \in \mathcal{K} \}``. After solving the scaled problem the original solution is obtained by reversing the scaling:
```math
   x = D\hat{x}, \quad s = E^{-1}\hat{s}, \quad y = E\hat{y}.
```
To obtain the scaling matrices ``D`` and ``E`` we use a modified Ruiz equilibration algorithm which involves a certain number of scaling iterations to equilibrate the column norms of the data matrices ``P`` and ``A``. The number of these iterations can be adjusted by the user with `scaling`. To disable the scaling step set `scaling = 0`.

## Termination criteria
The COSMO algorithm can terminate for four reasons:
* The maximum number of allowed iterations has been reached. The user can specify this value in the solver settings with `max_iter`.
* The solver runtime reaches the time limit specified by the user (`time_limit`).
* COSMO detects an infeasible problem.
* The iterates fulfil the termination criteria for convergence.

COSMO uses the primal and dual residuals of the problem to determine if the algorithm has converged. The primal and dual residuals are given by:
```math
\begin{align*}
r_p &= Ax + s -b,\\
r_d &= Px + q + A^\top y.
\end{align*}
```
The solver terminates when the ``\infty``-norms of the residuals lie below a specified tolerance. COSMO uses the sum of an absolute and relative tolerance term:
```math
\begin{align*}
  \|r_p^k \|_\infty &\leq \epsilon_{\mathrm{abs}} + \epsilon_{\mathrm{rel}} \, \text{max} \left\{ \|Ax^k \|_\infty,\|s^k\|_\infty, \|b\|_\infty \right\},\\
   \|r_d^k\|_\infty &\leq \epsilon_{\mathrm{abs}} + \epsilon_{\mathrm{rel}} \, \text{max} \left\{\|Px^k\|_\infty,\|q\|_\infty, \|A^\top y^k\|_\infty \right\}.
\end{align*}
```
The absolute and relative tolerances ``\epsilon_{\mathrm{abs}}``and ``\epsilon_{\mathrm{rel}}`` can be set by the user by specifying `eps_abs` and `eps_rel`. Furthermore, the user can adjust the number of iterations after which the convergence criteria are checked (`check_termination`).

## Infeasibility detection
COSMO uses conditions based on separating hyperplanes to detect infeasible problems. The conditions for COSMO's problem format have been developed in [1]. Define the convex set ``\mathcal{C} = \mathcal{-K} + \{b\}`` then we can use the following infeasibility conditions:
```math
\begin{align*}
\mathcal{P} &= \left\{x \in \mathbb{R}^n \mid  Px = 0, \, Ax \in \mathcal{C}^{\infty}, \, \langle  q,x \rangle < 0  \right\} \\
\mathcal{D} &= \left\{y \in \mathbb{R}^m \mid  A^\top  y  = 0,  \, \sigma_\mathcal{C}(y) < 0 \right\}.
\end{align*}
```
The existence of some ``y \in \mathcal{D}`` is a certificate that the problem is primal infeasible, while the existence of some ``x \in \mathcal{P}`` is a certificate for dual infeasibility. COSMO regularly checks above conditions to detect infeasible problems. If the detection is successful, the solver terminates and returns the status codes `:Primal_infeasible` or `:Dual_infeasible`. COSMO checks the conditions every `check_infeasibility` iterations, which can be adjusted by the user.

### References

[1] Banjac, G. et al. *Infeasibility detection in the alternating direction method of multipliers for convex optimization*. [Preprint](http://people.ee.ethz.ch/~gbanjac/pdfs/admm_infeas.pdf), 2017.


