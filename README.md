# OSQP Solver - Pure Julia Implementation
This is a pure Julia implementation of the OSQP solver available on [**Github**](https://github.com/oxfordcontrol/osqp). For more information read the documentation [**osqp.readthedocs.io**](http://osqp.readthedocs.io/).

The solver can handle problems of the form
```julia
minimize        1/2 x' P x + q' x

subject to      l <= A x <= u
```
where `x in R^n` is the optimization variable. The objective function is defined by a positive semidefinite matrix `P in S^n_+` and vector `q in R^n`. The linear constraints are defined by matrix `A in R^{m x n}` and vectors `l in R^m U {-inf}^m`, `u in R^m U {+inf}^m`.

## Usage
The solver can be used by including  **osqp_julia.jl** into your project. A short example is provided below:

```julia
# include solver file
include("osqp_julia.jl")
using OSQPSolver

# define sample problem
P = [4.0 1.0; 1.0 2.0]
q = [1.0;1.0]
A = [1.0 1.0; 1.0 0.0; 0.0 1.0]
l = [1.0;0.0;0.0]
u = [1.0;0.7;0.7]
settings = qpSettings(rho=1.0,verbose=true)

# solve QP problem
res = solveOSQP(P,q,A,l,u,settings)
```
## Contact
Send an email to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk)
