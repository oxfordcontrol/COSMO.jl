# Quadratic Objective Conic Solver (QOCS) - Pure Julia Implementation
This is a pure Julia implementation of the QOCS solver. It solves convex optimization problems of the following form:
```
min trace(X,PX)+trace(Q,X) 
s.t. trace(A_i*X) = b_i, i=1,...,m
     X in K
```
where K is a composite cone. The problem data has to be transformed into the following vector matrix description of the problem:
```
min 1/2 x'Px + q'x 
s.t. Ax = b, x in S+
```
## Installation / Usage
- Solver was written for Julia v0.6
- Clone repository to local machine
- Include the `Solver.jl` and `Projections.jl` files into your project and load the `OSSDP module`.
- Consider the following example:

```julia
include("../Projections.jl")
include("../Solver.jl")
using OSSDP


# Problem Data
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

# Reformat data to fit vector matrix input style
P = zeros(9,9)
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]

# set solver parameters
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true)

# solve problem
result,dbg = solveSDP(P,q,A,b,settings)

# Check results
@show(reshape(result.x,3,3))
@show(result.cost)
```
## Settings
Settings can be specified using the sdpSettings struct. The following settings are available:

Argument | Description | Values (default)
--- | --- | ---
rho | ADMM rho step | 1.0
sigma | ADMM sigma step | 10.0
alpha | Relaxation parameter | 1.6
eps_abs | Absolute residual tolerance | 1e-3
eps_rel | Relative residual tolerance | 1e-3
eps_prim_inf | Primal infeasibility tolerance | 1e-4
eps_dual_inf | Dual infeasibility tolerance | 1e-4
max_iter | Maximum number of iterations | 2500
verbose | Verbose printing | false
checkTermination | Check termination interval | 1
scaling | Number of scaling iterations | 10

## Tasks / Future Work
The current tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues):exclamation:

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!	
