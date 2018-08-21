![QOCS Logo](https://github.com/migarstka/QOCS_assets/blob/master/QOCS_logo.png)
This repository hosts a Julia implementation of the QOCS solver. It solves convex optimization problems of the following form:
```
min 1/2 x'Px + q'x
s.t. Ax + s = b, s in C
```
with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex set C is a composition of convex sets and cones. By default QOCS supports the zero cone, the non-negative orthant, second order cones and positive semidefinite cones. Further convex sets can be added by the user.

## Installation / Usage
- The Solver was written for Julia v0.7/1.0
- Clone package to local machine: `Pkg.clone("https://github.com/oxfordcontrol/QOCS.jl")`
- Load the package with `using QOCS`
- Consider the following example:

### Example
```julia

using QOCS, Test, LinearAlgebra

# Linear program example
# min c'x
# s.t. Ax <= b
#      x >= 1,  x2 =>5, x1+x3 => 4
c = [1; 2; 3; 4]
A = Matrix(1.0I,4,4)
b = [10; 10; 10; 10]

# create augmented matrices
Aa = -[A;-Matrix(1.0I,4,4);0 -1 0 0;-1 0 -1 0]
ba = vec([b; -ones(4,1);-5;-4])
P = zeros(size(A,2),size(A,2))

constraint1 = QOCS.Constraint(Aa,ba,QOCS.Nonnegatives())

# define example problem
settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,check_termination=1,eps_abs = 1e-6, eps_rel = 1e-6)

model = QOCS.Model()
assemble!(model,P,c,[constraint1])

res, = QOCS.optimize!(model,settings);

@testset "Linear Problem" begin
  @test isapprox(res.x[1:4],[3;5;1;1], atol=1e-2, norm=(x -> norm(x,Inf)))
  @test isapprox(res.cost,20.0, atol=1e-2)
end
nothing
```


## Settings
Settings can be specified using the `QOCS.Settings` struct. The following settings are available:

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
verboseTiming | Verbose timing | false
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
timelimit | set solver time limit in s | 0

For more low-level settings, see the Settings definition in `/src/Types.jl`.

## Result
After attempting to solve the problem, QOCS will return a result object with the following fields:

Fieldname | Type | Description
---  | --- | ---
x | Vector{Float64}| Primal variable
y | Vector{Float64}| Dual variable
s | Vector{Float64}| (Primal) set variable
objVal | Float64 | Objective value
iter | Int64 | Number of iterations
status | Symbol | Solution status
info | QOCS.ResultInfo | Struct with more information
times | QOCS.ResultTimes | Struct with several measured times

### Status Codes
QOCS will return one of the following statuses:

Status Code  | Description
---  | ---
:Solved | A optimal solution was found
:Unsolved | Default value
:Max_iter_reached | Solver reached iteration limit (set with `Settings.max_iter`)
:Time_limit_reached | Solver reached time limit (set with `Settings.timelimit`)
:Primal_infeasible | Problem is primal infeasible
:Dual_infeasible | Problem is dual infeasible


### Timings
If `settings.verboseTiming` is set to `true` QOCS will report the following times in `result.times`:

Time Name  | Description
---  | ---
solverTime | Total time used to solve the problem
setupTime | Setup time = graphTime + factorTime
graphTime | Time used to perform chordal decomposition
factorTime | Time used to vector the system of linear equations
iterTime | Time spent in iteration loop
projTime | Time spent in projection functions
postTime | Time used for post processing

It holds: `solverTime` = `setupTime`+ `iterTime`, `setupTime` = `graphTime`+ `factorTime`, `projTime` subset of `iterTime`.


## Test problems
A set of benchmark problems with conic constraints has been collected and made available here:
[https://github.com/migarstka/SDP_Benchmark_Problems](https://github.com/migarstka/SDP_Benchmark_Problems)

## Tasks / Future Work
The current tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!
