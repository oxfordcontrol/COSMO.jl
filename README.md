![QOCS Logo](https://github.com/migarstka/QOCS_assets/blob/master/QOCS_logo.png)
This repository hosts a Julia implementation of the QOCS solver. It solves convex optimization problems of the following form:
```
min 1/2 x'Px + q'x
s.t. Ax + s = b, s in K
```
with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex cone K is a composition of the zero cone, the non-negative orthant, a set of second order cones, and a set of positive semidefinite cones. The dimension of the cones have to be specified using the `Cone` type (`K.f::Int`: number of zero cone variables, `K.l::Int`: number of nonnegative components, `K.s::Array{Int}`: number of variables in each second-order cone, `K.q::Array{Int}`: number of variables in each psd cone).

## Installation / Usage
- The Solver was written for Julia v0.6
- Clone repository to local machine
- Include `../src/QOCS.jl` into your project and load the `QOCS` module.
- Consider the following example:

```julia
workspace()
include("../src/QOCS.jl")

using Base.Test
using QOCS

# Linear program example
# min c'x
# s.t. Ax <= b
#      x >= 1,  x2 =>5, x1+x3 => 4
c = [1; 2; 3; 4]
A = eye(4)
b = [10; 10; 10; 10]

# create augmented matrices
Aa = [A;-eye(4);0 -1 0 0;-1 0 -1 0]
ba = [b; -ones(4,1);-5;-4]
P = zeros(size(A,2),size(A,2))

# define cone dimensions
K = QOCS.Cone(0,10,[],[])

# adjust solver settings
settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,check_termination=1,scaling = 0,eps_abs = 1e-6, eps_rel = 1e-6)

# solve problem
res,ws  = QOCS.solve(P,c,Aa,ba,K,settings);

# test against known solution
@testset "Linear Problem" begin
  @test isapprox(res.x[1:4],[3;5;1;1], atol=1e-2, norm=(x -> norm(x,Inf)))
  @test isapprox(res.cost,20.0, atol=1e-2)
end
```
## Test problems
A set of benchmark problems with conic constraints has been collected and made available here:
[https://github.com/migarstka/SDP_Benchmark_Problems](https://github.com/migarstka/SDP_Benchmark_Problems)
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
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
timelimit | set solver time limit in s | 0

For more low-level settings, see the Settings definition in `/src/Types.jl`.

## Status Codes
After attempting to solve the problem, QOCS will return one of the following statuses:

Status Code  | Description
---  | ---
:Solved | A optimal solution was found
:Unsolved | Default value
:Max_iter_reached | Solver reached iteration limit (set with Settings.max_iter)
:Time_limit_reached | Solver reached time limit (set with Settings.timelimit)
:Primal_infeasible | Problem is primal infeasible
:Dual_infeasible | Problem is dual infeasible


## Tasks / Future Work
The current tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!
