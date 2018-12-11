![COSMO Logo](https://github.com/migarstka/COSMO_assets/blob/master/COSMO_logo.png)
This repository hosts a Julia implementation of the COSMO solver. It solves convex optimization problems of the following form:
```
min 1/2 x'Px + q'x
s.t. Ax + s = b, s in C
```
with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex set `C` is a composition of convex sets and cones. By default COSMO supports the zero cone, the non-negative orthant, second order cones and positive semidefinite cones. Further convex sets can be added by the user.

## Installation / Usage
- The Solver was written for Julia v1.0
- Add package via the package manager (type `]`): `add https://github.com/oxfordcontrol/COSMO.jl`
- Make the package available with `using COSMO`
- Consider the following example:

### Example
```julia

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
c1 = COSMO.Constraint(A, b, COSMO.Nonnegatives)
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


## Settings
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
time_limit | set solver time limit in s | 0

For more low-level settings, see the Settings definition in `/src/Types.jl`.

## Result
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
If `settings.verbose_timing` is set to `true` COSMO will report the following times in `result.times`:

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
`proj_time` subset of `iter_time`.


## Test problems
A set of benchmark problems with conic constraints has been collected and made available here:
[https://github.com/migarstka/SDP_Benchmark_Problems](https://github.com/migarstka/SDP_Benchmark_Problems)

## Tasks / Future Work
The current tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!
