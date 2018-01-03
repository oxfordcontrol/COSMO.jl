workspace()
include("../osqp_julia.jl")
using OSQPSolver

# define example problem
P = [4.0 1.0; 1.0 2.0]
q = [1.0;1.0]
A = [1.0 1.0; 1.0 0.0; 0.0 1.0]
l = [1.0;0.0;0.0]
u = [1.0;0.7;0.7]
settings = qpSettings(rho=1.0,verbose=true)

# solve QP problem
res = solveOSQP(P,q,A,l,u,settings)
nothing;
