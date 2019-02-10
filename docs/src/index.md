__COSMO.jl__ is a Julia implementation of the _Conic Operator Splitting Method_. The underlying ADMM-algortihm is well-suited for large convex conic problems. COSMO solves the following problem:

```math
\begin{array}{ll} \mbox{minimize} & \textstyle{\frac{1}{2}}x^\top Px + q^\top x\\ \mbox{subject to} & Ax + s = b \\ & s \in \mathcal{C}, \end{array}
```

with decision variables ``x \in \mathrm{R}^n``, ``s \in \mathrm{R}^m`` and data matrices ``P=P^\top \succeq 0``, ``q \in \mathrm{R}^n``, ``A \in \mathrm{R}^{m \times n}``, and ``b \in \mathrm{R}^m``. The convex set ``\mathcal{C}``
 is a composition of convex sets and cones.

## Features

* COSMO solves linear programs, quadratic programs, second-order cone programs and semidefinite programs
* Semi-definite programs with quadratic objective functions are natively supported
* Infeasible problems are detected without a homogeneous self-dual embedding of the problem
* Support for MathOptInterface and the upcoming JuMP `v0.19` release, which allows you to describe your problem in JuMP.
* COSMO tries to decompose large structured PSD constraints using chordal decomposition techniques. This often results in a significant speedup compared to the original problem. (_this feature is not yet available in Julia `v1.0`_)
* COSMO supports warm-starting of the decision variables.

## Installation
COSMO can be installed using the Julia package manager for Julia `v1.0` and higher. Inside the Julia REPL, type `]` to enter the Pkg REPL mode then run

`pkg> add COSMO`

If you want to install the latest version from master run

`pkg> add COSMO#master`

## Quick Example

## Credits

The following people are involved in the development of COSMO:
* Michael Garstka (main development)
* Nikitas Rontsis (algorithm performance)
* Paul Goulart (code architecture, maths and algorithms)
* Mark Cannon (maths and algorithms)
\*all contributors are affiliated with the University of Oxford.

If this project is useful for your work please consider
* [Citing](citing.md) the relevant paper
* Leaving a star on the [GitHub repository](https://github.com/oxfordcontrol/COSMO.jl)





## Licence
COSMO.jl is licensed under the Apache License 2.0. For more details click [here](https://github.com/oxfordcontrol/COSMO.jl/blob/master/LICENSE.md).