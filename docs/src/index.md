__COSMO.jl__ is a Julia implementation of the _Conic Operator Splitting Method_. The underlying ADMM-algorithm is well-suited for large convex conic problems. COSMO solves the following problem:

```math
\begin{array}{ll} \mbox{minimize} & \textstyle{\frac{1}{2}}x^\top Px + q^\top x\\ \mbox{subject to} & Ax + s = b \\ & s \in \mathcal{K}, \end{array}
```

with decision variables ``x \in \mathbb{R}^n``, ``s \in \mathbb{R}^m`` and data matrices ``P=P^\top \succeq 0``, ``q \in \mathbb{R}^n``, ``A \in \mathbb{R}^{m \times n}``, and ``b \in \mathbb{R}^m``. The convex set ``\mathcal{K}``
 is a composition of convex sets and cones.

## Features

* __Versatile__: COSMO solves linear programs, quadratic programs, second-order cone programs, semidefinite programs and problems involving exponential and power cones
* __Quad SDPs__: Positive semidefinite programs with quadratic objective functions are natively supported
* __Infeasibility detection__: Infeasible problems are detected without a homogeneous self-dual embedding of the problem
* __JuMP support__: COSMO supports MathOptInterface and JuMP, which allows you to model your problem in JuMP
* __Chordal decomposition__: COSMO tries to decompose large structured PSD constraints using chordal decomposition techniques. This often results in a significant speedup compared to the original problem. (_this feature is not active by default_)
* __Warm starting__: COSMO supports warm starting of the decision variables
* __Open Source__: Our code is available on [GitHub](https://github.com/oxfordcontrol/COSMO.jl) and distributed under the Apache 2.0 Licence

## Installation
COSMO can be installed using the Julia package manager for Julia `v1.0` and higher. Inside the Julia REPL, type `]` to enter the Pkg REPL mode then run

`pkg> add COSMO`

If you want to install the latest version from master run

`pkg> add COSMO#master`

## Quick Example
Consider the following 2x2 semidefinite program with decision variable `X`:
```math
\begin{array}{ll} \mbox{minimize} &  \text{tr}(CX)\\
\mbox{subject to} &  \text{tr}(A X) = b \\
                  &  X \succeq 0,
\end{array}
```
with problem data `A`, `b` and `C`:
```math
A = \begin{bmatrix} 1 & 0 \\ 5 & 2\end{bmatrix},
C = \begin{bmatrix} 1 & 2 \\ 0 & 2\end{bmatrix},
b = 4.
```
where `tr` denotes the trace of a matrix.
We can solve this problem either using COSMO's interface:
```julia
using COSMO, LinearAlgebra

C =  [1. 2; 0 2]
A = [1. 0; 5 2]
b = 4;

model = COSMO.Model();

# define the cost function
P = zeros(4, 4)
q = vec(C)

# define the constraints
# A x = b
cs1 = COSMO.Constraint(vec(A)', -b, COSMO.ZeroSet)
# X in PSD cone
cs2 = COSMO.Constraint(Matrix(1.0I, 4, 4), zeros(4), COSMO.PsdCone)
constraints = [cs1; cs2]

# assemble and solve the model
assemble!(model, P, q, constraints)
result = COSMO.optimize!(model);

X_sol = reshape(result.x, 2, 2)
obj_value = result.obj_val
```

or we can describe the problem using `JuMP` and use COSMO as the backend solver:
```julia
using COSMO, JuMP, LinearAlgebra

C =  [1 2; 0 2]
A = [1 0; 5 2]
b = 4;

m = Model(with_optimizer(COSMO.Optimizer));
@variable(m, X[1:2, 1:2], PSD)
@objective(m, Min, tr(C * X));
@constraint(m, tr(A * X) == b);
JuMP.optimize!(m);

status = JuMP.termination_status(m)
X_sol = JuMP.value.(X)
obj_value = JuMP.objective_value(m)
```

## Presentation Video
A video of the presentation at JuMP-dev is available here:
```@raw html
<a href="http://www.youtube.com/watch?feature=player_embedded&v=VTxPqEzzqlU" target="_blank"><img src="https://raw.githubusercontent.com/migarstka/COSMO_assets/master/jump_dev_video.png" alt="JuMP dev presentation video" width="400" border="10"/></a>
```

## Credits

The following people are involved in the development of COSMO:
* [Michael Garstka](https://migarstka.github.io) (main development)
* [Nikitas Rontsis](https://github.com/nrontsis) (algorithm performance)
* [Paul Goulart](http://users.ox.ac.uk/~engs1373/) (code architecture, maths and algorithms)
* [Mark Cannon](https://markcannon.github.io) (maths and algorithms)
\*all contributors are affiliated with the [University of Oxford](http://www2.eng.ox.ac.uk/control).

If this project is useful for your work please consider
* [Citing](citing.md) the relevant paper
* Leaving a star on the [GitHub repository](https://github.com/oxfordcontrol/COSMO.jl)


## Licence
COSMO.jl is licensed under the Apache License 2.0. For more details click [here](https://github.com/oxfordcontrol/COSMO.jl/blob/master/LICENSE.md).
