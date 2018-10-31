workspace()
include("../../../src/Solver.jl")
using JLD
using Base.Test
using OSSDP
using JuMP, Mosek

problemName = "truss1.jld"
data = JLD.load("../sdplib/"*problemName)
F = data["F"]
c = data["c"]
m = data["m"]
n = data["n"]
optVal = data["optVal"]

# Rewrite problem to OSSDP compatible format:
# min   1/2 x'Px + q'x
# s.t.  Ax + s = b
#       s in K
# -----------------------
# primal form
# -----------------------
P = zeros(m,m)
q = c
A = zeros(n^2,m)
for iii = 1:m
  A[:,iii] = -vec(F[iii+1])
end
b = -vec(F[1])
Kf = 0
Kl = 0
Kq = []
Ks = [n^2]


# solve with COSMO
K = Cone(Kf,Kl,Kq,Ks)
settings = OSSDPSettings(max_iter=3000,verbose=false,check_termination=1,checkInfeasibility=50,scaling = 10 ,scaleFunc=2,adaptive_rho=true,eps_abs=1e-6,eps_rel=1e-6)
res,nothing = OSSDP.solve(P,q,A,b,K,settings);

# solve with MOSEK
model = Model(solver=MosekSolver())
@variable(model, x[1:m])
@variable(model, S[1:n,1:n],SDP)
@objective(model, Min, q'*x)
s = vec(S)
@constraint(model, A*x+s .== b)
status = JuMP.solve(model)

@testset "Test SDPLib sample problem" begin
    @test abs( res.cost - getobjectivevalue(model) ) < 1e-2
end



