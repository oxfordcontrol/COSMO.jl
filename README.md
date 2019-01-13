
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
  <a href="https://opensource.org/licenses/Apache-2.0"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg"></a>
</p>



<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#examples">Examples</a> •
  <a href="#interface">Interface</a> •
  <a href="NEWS.md">News</a> •
  <a href="#contributing">Contributing</a> •
  <a href="#contact">Contact</a>
</p>

This is a Julia implementation of the _Conic operator splitting method_ (COSMO) solver. It can solve large convex conic optimization problems of the following form:
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\begin{array}{ll}&space;\mbox{minimize}&space;&&space;\textstyle{\frac{1}{2}}x^\top&space;Px&space;&plus;&space;q^\top&space;x\\&space;\mbox{subject&space;to}&space;&&space;Ax&space;&plus;&space;s&space;=&space;b&space;\\&space;&&space;s&space;\in&space;\mathcal{C},&space;\end{array}" title="\begin{array}{ll} \mbox{minimize} & \textstyle{\frac{1}{2}}x^\top Px + q^\top x\\ \mbox{subject to} & Ax + s = b \\ & s \in \mathcal{C}, \end{array}"/>
</p>

with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex set `C` is a composition of convex sets and cones.

## Features
By default COSMO supports the zero cone, the non-negative orthant, second order cones and positive semidefinite cones. COSMO allows you to:
- solve semidefinite programs with quadratic objective functions directly
- detect infeasible problems without a homogeneous self-dual embedding of the problem
- describe your optimisation problem using [JuMP](https://github.com/JuliaOpt/JuMP.jl) (COSMO requires JuMP v0.19-beta)
- use chordal decomposition techniques to decompose chordally structured SDPs
- define your own convex sets for constraints



## Installation
- The solver is written for Julia `v1.0`
- Add the package via the package manager (type `]`): `>> add COSMO`
- Make the package available in your project with `>> using COSMO`

## Examples

### Using JuMP
We consider the problem of finding the closest correlation matrix X, i.e. PSD and ones on the diagonal, to a random matrix C.
```julia
using COSMO, JuMP, LinearAlgebra, SparseArrays, Test, Random
rng = Random.MersenneTwister(12345);

# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

# create a random test matrix C
n = 8
C = -1 .+ rand(rng, n, n) .* 2;
c = vec(C);

# define problem in JuMP
q = -vec(C);
r = 0.5 * vec(C)' * vec(C);
m = Model(with_optimizer(COSMO.Optimizer, verbose=true, eps_abs = 1e-4));
@variable(m, X[1:n, 1:n], PSD);
x = vec(X);
@objective(m, Min, 0.5 * x' * x  + q' * x + r)
for i = 1:n
  @constraint(m, X[i, i] == 1.)
end

# solve and get results
status = JuMP.optimize!(m)
obj_val = JuMP.objective_value(m)
X_sol = JuMP.value.(X)
```

### Using the native solver interface
```julia
using COSMO, LinearAlgebra, SparseArrays, Test

# Linear program example
# min c'x
# s.t. Ax <= b
#      x >= 1,  x2 >= 5, x1+x3 >= 4
# with data matrices
c = [1; 2; 3; 4.]
A = Matrix(1.0I, 4, 4)
b = [10; 10; 10; 10]
n = 4
# -------------------
# create constraints A * x + b in set
# -------------------
# Ax <= b
c1 = COSMO.Constraint(-A, b, COSMO.Nonnegatives)
# x >= 1
c2 = COSMO.Constraint(Matrix(1.0I, n, n), -ones(n), COSMO.Nonnegatives)
# x2 >= 5
c3 = COSMO.Constraint(1, -5, COSMO.Nonnegatives, n, 2:2)
# x1 + x3 >= 4
c4 = COSMO.Constraint([1 0 1 0], -4, COSMO.Nonnegatives)

# -------------------
# define cost function
# -------------------
P = spzeros(4, 4)
q = c

# -------------------
# assemble solver model
# -------------------
settings = COSMO.Settings(max_iter=2500, verbose=true, eps_abs = 1e-4, eps_rel = 1e-5)
model = COSMO.Model()
assemble!(model, P, q, [c1; c2; c3; c4], settings)
res = COSMO.optimize!(model);

@testset "Linear Problem" begin
  @test isapprox(res.x[1:4], [3; 5; 1; 1], atol=1e-2, norm = (x -> norm(x, Inf)))
  @test isapprox(res.obj_val, 20.0, atol=1e-2)
end
```


### Test problems
A set of benchmark problems with conic constraints is available here:
[https://github.com/migarstka/SDP_Benchmark_Problems](https://github.com/migarstka/SDP_Benchmark_Problems)

## Interface

### Settings
Settings can be specified using the `COSMO.Settings` struct. The following settings are available:

Argument | Description | Values (default)
--- | --- | ---
rho | ADMM rho step | 0.1
sigma | ADMM sigma step | 1e-6.
alpha | Relaxation parameter | 1.6
eps_abs | Absolute residual tolerance | 1e-4
eps_rel | Relative residual tolerance | 1e-4
eps_prim_inf | Primal infeasibility tolerance | 1e-4
eps_dual_inf | Dual infeasibility tolerance | 1e-4
max_iter | Maximum number of iterations | 2500
verbose | Verbose printing | false
verbose_timing | Verbose timing | false
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
time_limit | set solver time limit in s | 0.0

For more low-level settings, see the `Settings` type definition in `/src/types.jl`.

### Result
After attempting to solve the problem, COSMO will return a result object with the following fields:

Fieldname | Type | Description
---  | --- | ---
x | Vector{Float64}| Primal variable
y | Vector{Float64}| Dual variable
s | Vector{Float64}| (Primal) set variable
obj_val | Float64 | Objective value
iter | Int64 | Number of iterations
status | Symbol | Solution status
info | COSMO.ResultInfo | Struct with more information
times | COSMO.ResultTimes | Struct with several measured times

### Status Codes
COSMO will return one of the following statuses:

Status Code  | Description
---  | ---
:Solved | A optimal solution was found
:Unsolved | Default value
:Max_iter_reached | Solver reached iteration limit (set with `Settings.max_iter`)
:Time_limit_reached | Solver reached time limit (set with `Settings.timelimit`)
:Primal_infeasible | Problem is primal infeasible
:Dual_infeasible | Problem is dual infeasible


### Timings
If `settings.verbose_timing` is set to `true`, COSMO will report the following times in `result.times`:

Time Name  | Description
---  | ---
solver_time | Total time used to solve the problem
setup_time | Setup time = graphTime + factorTime
graph_time | Time used to perform chordal decomposition
factor_time | Time used to factor the system of linear equations
iter_time | Time spent in iteration loop
proj_time | Time spent in projection functions
post_time | Time used for post processing

It holds:
`solver_time` = `setup_time`+ `iter_time` + `post_time`,
`setup_time` = `graph_time`+ `factor_time`,
`proj_time` is a subset of `iter_time`.

## Contributing
- Current issues, tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:. Please report any issues or bugs that you encounter.
- Contributions are always welcome. Our style guide can be found [here](https://github.com/oxfordcontrol/COSMO.jl/wiki/Code-Style-Guide).
- As an open source project we are also interested in any projects and applications that use COSMO. Please let us know! 

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!
