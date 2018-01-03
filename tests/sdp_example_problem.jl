# Matrix minimization problem
workspace()
using JuMP,SCS


A = [1.0 0.0;0.0 0.0; 0.0 1.0;0 -1;0 0;1 0; 0 1; 0 0; 0 1;0 2;1 0; 0 0;0 -1;0 0;0 0;1 0]
b = [0.0;0;-2;0;0;0;0;-2;-2;0;0;0;0;-2;0;0]
c = [1;0]
m = size(A,1)
n = size(A,2)

mod = Model(solver=SCSSolver())
@variable(mod,x)
@variable(mod,s)
@variable(mod,z[1:16])
@objective(mod,Min,s)
@constraint(mod,A*[s;x]-z .== b)
@SDconstraint(m, X >= eye(n))

status = solve(mod)