
<div style="display: flex; justify-content: flex-end" margin=0px>
    <img src="https://github.com/migarstka/COSMO_assets/blob/master/star_badge_2.png" align="right" width=6%>
</div>
<h1 align="center" margin=0px>
  <img src="https://github.com/migarstka/COSMO_assets/blob/master/cosmo_rocket_with_convex_set.png" width=70%>
</h1>
<p align="center">
  <a href="https://travis-ci.org/oxfordcontrol/COSMO.jl"><img src="https://travis-ci.org/oxfordcontrol/COSMO.jl.svg?branch=master"></a>
  <a href="https://codecov.io/gh/oxfordcontrol/COSMO.jl"><img src="https://codecov.io/gh/oxfordcontrol/COSMO.jl/branch/master/graph/badge.svg"></a>
  <a href="https://oxfordcontrol.github.io/COSMO.jl/stable"><img src="https://img.shields.io/badge/Documentation-stable-purple.svg"></a>
  <a href="https://opensource.org/licenses/Apache-2.0"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg"></a>
  <a href="https://github.com/oxfordcontrol/COSMO.jl/releases"><img src="https://img.shields.io/badge/Release-v0.7.8-blue.svg"></a>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="NEWS.md">News</a> •
  <a href="#citing-">Citing</a> •
  <a href="#contributing">Contributing</a>
</p>

This is a Julia implementation of the _Conic operator splitting method_ (COSMO) solver. It can solve large convex conic optimization problems of the following form:
<p align="center">
<img src="https://github.com/migarstka/COSMO_assets/blob/master/cosmo_format.png" width=220px>
</p>


with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex set `K` is a composition of convex sets and cones.

__For more information take a look at the COSMO.jl Documentation ([stable](https://oxfordcontrol.github.io/COSMO.jl/stable) |  [dev](https://oxfordcontrol.github.io/COSMO.jl/dev)).__

## Features

* __Versatile__: COSMO solves linear programs, quadratic programs, second-order cone programs, semidefinite programs and problems involving exponential and power cones
* __Quad SDPs__: Positive semidefinite programs with quadratic objective functions are natively supported
* __Infeasibility detection__: Infeasible problems are detected without a homogeneous self-dual embedding of the problem
* __JuMP / Convex.jl support__: We provide an interface to MathOptInterface (MOI), which allows you to describe your problem in [JuMP](https://github.com/JuliaOpt/JuMP.jl) and [Convex.jl](https://github.com/JuliaOpt/Convex.jl).
* __Warm starting__: COSMO supports warm starting of the decision variables
* __Custom sets and linear solver__: Customize COSMO's components by defining your own convex constraint sets and by choosing from a number of direct and indirect linear system solvers, e.g. [QDLDL](https://github.com/oxfordcontrol/qdldl), [Pardiso](https://github.com/JuliaSparse/Pardiso.jl), [Conjugate Gradient](https://juliamath.github.io/IterativeSolvers.jl/dev/) and [MINRES](https://juliamath.github.io/IterativeSolvers.jl/dev/)
* __Arbitrary precision types__: You can solve problems with any floating point precision.
* __Open Source__: Our code is free to use and distributed under the Apache 2.0 Licence
* __Chordal decomposition__: COSMO tries to decompose large structured PSD constraints using chordal decomposition techniques. This often results in a significant speedup compared to the original problem.
* __Smart clique merging__: After an initial decomposition of a structured SDP, COSMO recombines overlapping cliques/blocks to speed up the algorithm.
<div align="center" margin=0px>
  <img src="docs/src/assets/example_decomposition.gif" width=45%>
</div>

## Installation
- `COSMO` can be added via the Julia package manager (type `]`): `pkg> add COSMO`

## Citing
If you find COSMO useful in your project, we kindly request that you cite the following paper:
```
@InProceedings{garstka_2019,
  author        = {Michael Garstka and Mark Cannon and Paul Goulart},
  title         = {{COSMO}: A conic operator splitting method for large convex problems},
  booktitle     = {European Control Conference},
  year          = {2019},
  location      = {Naples, Italy},
  doi            = {10.23919/ECC.2019.8796161},
  eprint        = {1901.10887},
  url           = {https://arxiv.org/abs/1901.10887},
  archiveprefix = {arXiv},
  keywords      = {Mathematics - Optimization and Control},
  primaryclass  = {math.OC},
}
```
A preprint can be downloaded [here](https://arxiv.org/abs/1901.10887).



## Contributing
- Contributions are always welcome. Our style guide can be found [here](https://github.com/oxfordcontrol/COSMO.jl/wiki/Code-Style-Guide).
- Current issues, tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues). Please report any issues or bugs that you encounter.
- As an open source project we are also interested in any projects and applications that use COSMO. Please let us know!

## Python - Interface
COSMO can also be called from Python. Take a look at: [cosmo-python](https://github.com/oxfordcontrol/cosmo-python)

## Licence 🔍
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.
