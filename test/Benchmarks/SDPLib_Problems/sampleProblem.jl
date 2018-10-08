include("../../../src/QOCS.jl")
using JLD2
using Test
using Main.QOCS
using SparseArrays
using JuMP, SCS

filename = "../sdplib/truss1.jld"
# filename = "../sdplib/arch0.jld"
@load filename m n nblocks blockvec c F optVal

# Rewrite problem to OSSDP compatible format:
# min   1/2 x'Px + q'x
# s.t.  Ax + s = b
#       s in K
# -----------------------
# primal form
# -----------------------
# Problem data
P = zeros(m,m)
q = c
A = zeros(n^2,m)
for iii = 1:m
  A[:,iii] = -vec(F[iii+1])
end
b = -vec(F[1])
######

# solve with QOCS
constraint1 = QOCS.Constraint(A,b,QOCS.PositiveSemidefiniteCone())
settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=true,check_termination=1,eps_abs = 1e-6, eps_rel = 1e-6)
model = QOCS.Model()
assemble!(model,P,q,(constraint1))
res = QOCS.optimize!(model,settings);

# solve with MOSEK
mosek_model = JuMP.Model()
setsolver(mosek_model, SCSSolver())
@variable(mosek_model, x[1:m])
@variable(mosek_model, S[1:n,1:n],SDP)
@objective(mosek_model, Min, q'*x)
s = vec(S)
@constraint(mosek_model, A*x+s .== b)
status = JuMP.solve(mosek_model)

@testset "Test SDPLib sample problem" begin
    @test abs( res.objVal - getobjectivevalue(mosek_model) ) < 1e-2
end

nothing