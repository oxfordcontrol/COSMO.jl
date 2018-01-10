# OSSDP Solver - Pure Julia Implementation
This is a pure Julia implementation of the OSSDP solver. It solves semidefinite programs of the following form:
```
min trace(X,PX)+trace(Q,X) 
s.t. trace(A_i*X) = b_i, i=1,...,m
     X in S+
```
The problem data has to be transformed into the following vector matrix description of the problem:
```
min 1/2 x'Px + q'x 
s.t. Ax = b, x in S+
```

**Work in Progress!**
## Installation / Usage
- Clone repository to local machine
- include the ossdp.jl file into your project and load the OSSDP module.

```julia
include("../ossdp.jl")
using OSSDP


# Problem Data
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

# Reformat data to fit vector matrix input style
P = zeros(3,3)
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]

# set solver parameters
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true)

# solve problem
result,dbg = solveSDP(P,q,A,b,settings)

# Check results
@show(reshape(res.x,3,3))
@show(res.cost)
```

## Tasks / Future Work
The current tasks and future ideas are listed in [Issues](https://github.com/oxfordcontrol/ossdp/issues) :exclamation:

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!	
