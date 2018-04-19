# Quadratic Objective Conic Solver (QOCS) - Pure Julia Implementation
This is a pure Julia implementation of the QOCS solver. It solves convex optimization problems of the following form:
```
min_x 1/2 x'Px + q'x 
s.t. Ax + s = b, s in K
```
with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex cone K is a composition of the zero cone, the non-negative orthant, a set of second order cones, and a set of positive semidefinite cones. The dimension of the cones have to be specified using the `Cone` type (`K.f::Int`: number of zero cone variables, `K.l::Int`: number of nonnegative components, `K.s::Array{Int}`: number of variables in each second-order cone, `K.q::Array{Int}`: number of variables in each psd cone).

## Installation / Usage
- The Solver was written for Julia v0.6
- Clone repository to local machine
- Include `../src/Solver.jl` into your project and load the `OSSDP` and `OSSDPTypes` module.
- Consider the following example:

```julia
workspace()
include("../src/Solver.jl")

using Base.Test
using OSSDP, OSSDPTypes

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
K = OSSDPTypes.Cone(0,10,[],[])

# adjust solver settings
settings = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0,eps_abs = 1e-6, eps_rel = 1e-6)

# solve problem
res,ws  = OSSDP.solve(P,c,Aa,ba,K,settings);

# test against known solution
@testset "Linear Problem" begin
  @test isapprox(res.x[1:4],[3;5;1;1], atol=1e-2, norm=(x -> norm(x,Inf)))
  @test isapprox(res.cost,20.0, atol=1e-2)
end
```
## Settings
Settings can be specified using the `OSSDPSettings` struct. The following settings are available:

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
checkTermination | Check termination interval | 1
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | false
adaptive_rho_interval | Number of iterations after which rho is adapted | 40

For more low-level settings, see the OSSDPSettings type definition in `/src/Types.jl`.

## Tasks / Future Work
The current tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!	
