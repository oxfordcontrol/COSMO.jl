# User Guide
This user guide describes the basic structures and functions to define an optimisation problem, to solve the problem and to analyse the result. If you want to use `JuMP` to describe the problem, see the [JuMP Interface](@ref) section.

## Model
The problem data, user settings and workspace variables are all stored in a `Model`. To get started define an empty model:
```julia
m = COSMO.Model()
```
To initialize the model with an optimisation problem we need to define three more things:
* the objective function, i.e. the matrix `P` and the vector `q` in ``\frac{1}{2}x^\top P x + q^\top x``
* an array of constraints
* a `Settings` object that specifies how COSMO solves the problem _(optional)_

## Objective Function
To set the objective function of your optimisation problem simply define the square positive semidefinite matrix ``P \in \mathrm{R}^{n\times n} `` and the vector ``q \in \mathrm{R}^{n}``. You might have to transform your optimisation problem for this step.

## Constraints
```@docs
COSMO.Constraint
```

## Convex Sets

```@docs
COSMO.ZeroSet
COSMO.Nonnegatives
COSMO.SecondOrderCone
COSMO.PsdCone
COSMO.PsdConeTriangle
```

## Settings

Settings can be specified using the `COSMO.Settings` struct. The following settings are available:

Argument | Description | Values (default)
--- | --- | ---
rho | ADMM rho step | 0.1
sigma | ADMM sigma step | 1e-6
alpha | Relaxation parameter | 1.6
eps_abs | Absolute residual tolerance | 1e-4
eps_rel | Relative residual tolerance | 1e-4
eps\_prim\_inf | Primal infeasibility tolerance | 1e-4
eps\_dual\_inf | Dual infeasibility tolerance | 1e-4
max_iter | Maximum number of iterations | 2500
verbose | Verbose printing | false
verbose_timing | Verbose timing | false
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
time_limit | set solver time limit in s | 0.0

For more low-level settings, see the `Settings` type definition in `/src/settings.jl`.

## Results

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
:Max\_iter\_reached | Solver reached iteration limit (set with `Settings.max_iter`)
:Time\_limit\_reached | Solver reached time limit (set with `Settings.time_limit`)
:Primal\_infeasible | Problem is primal infeasible
:Dual\_infeasible | Problem is dual infeasible

### Timings
If `settings.verbose_timing` is set to `true`, COSMO will report the following times in `result.times`:

Time Name  | Description
---  | ---
solver_time | Total time used to solve the problem
setup_time | Setup time = graph\_time + factor\_time
graph_time | Time used to perform chordal decomposition
factor_time | Time used to factor the system of linear equations
iter_time | Time spent in iteration loop
proj_time | Time spent in projection functions
post_time | Time used for post processing

It holds:
`solver_time` = `setup_time`+ `iter_time` + `post_time`,
`setup_time` = `graph_time`+ `factor_time`,
`proj_time` is a subset of `iter_time`.


