
<div style="display: flex; justify-content: flex-end" margin=0px>
    <img src="https://github.com/migarstka/COSMO_assets/blob/master/star_badge_2.png" align="right" width=6%>
</div>
<h1 align="center" margin=0px>
  <br>
  <img src="https://github.com/migarstka/COSMO_assets/blob/master/COSMO_logo_only.png" width=40%>
  <br>
  <img src="https://github.com/migarstka/COSMO_assets/blob/master/COSMO_text_only.png" width=50%>
  <br>
</h1>
<p align="center">
  <a href="https://travis-ci.org/oxfordcontrol/COSMO.jl"><img src="https://travis-ci.org/oxfordcontrol/COSMO.jl.svg?branch=master"></a>
  <a href="https://codecov.io/gh/oxfordcontrol/COSMO.jl"><img src="https://codecov.io/gh/oxfordcontrol/COSMO.jl/branch/master/graph/badge.svg"></a>
  <a href="https://oxfordcontrol.github.io/COSMO.jl/stable"><img src="https://img.shields.io/badge/Documentation-stable-purple.svg"></a>
  <a href="https://opensource.org/licenses/Apache-2.0"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg"></a>
  <a href="https://github.com/oxfordcontrol/COSMO.jl/releases"><img src="https://img.shields.io/badge/Release-v0.5.0-blue.svg"></a>
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
<img src="https://latex.codecogs.com/gif.latex?\begin{array}{ll}&space;\mbox{minimize}&space;&&space;\textstyle{\frac{1}{2}}x^\top&space;Px&space;&plus;&space;q^\top&space;x\\&space;\mbox{subject&space;to}&space;&&space;Ax&space;&plus;&space;s&space;=&space;b&space;\\&space;&&space;s&space;\in&space;\mathcal{K},&space;\end{array}" title="\begin{array}{ll} \mbox{minimize} & \textstyle{\frac{1}{2}}x^\top Px + q^\top x\\ \mbox{subject to} & Ax + s = b \\ & s \in \mathcal{C}, \end{array}"/>
</p>

with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex set `K` is a composition of convex sets and cones.

__For more information check the [COSMO.jl Documentation](https://oxfordcontrol.github.io/COSMO.jl/stable).__

## Features

* __Versatile__: COSMO solves linear programs, quadratic programs, second-order cone programs, semidefinite programs and problems involving exponential and power cones
* __Quad SDPs__: Positive semidefinite programs with quadratic objective functions are natively supported
* __Infeasibility detection__: Infeasible problems are detected without a homogeneous self-dual embedding of the problem
* __JuMP support__: COSMO supports MathOptInterface and JuMP, which allows you to describe your problem in JuMP
* __Chordal decomposition__: COSMO tries to decompose large structured PSD constraints using chordal decomposition techniques. This often results in a significant speedup compared to the original problem.
* __Warm starting__: COSMO supports warm starting of the decision variables
* __Open Source__: Our code is free to use and distributed under the Apache 2.0 Licence

## Installation
- `COSMO` can be added via the Julia package manager (type `]`): `pkg> add COSMO`

## Citing 📃
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
- Current issues, tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:. Please report any issues or bugs that you encounter.
- As an open source project we are also interested in any projects and applications that use COSMO. Please let us know!

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.


