

workspace()
include("../../src/Solver.jl")
include("ConvertProblem.jl")
using Base.Test, OSSDP, JLD, MAT, OSQP

filename ="GENHS28"
data = matread("./MAT_FILES/$(filename).mat")

myP = zeros(2,2)

Q = data["Q"]
A = data["A"]
c = data["c"]
ru = data["ru"]
rl = data["rl"]
lb = data["lb"]
ub = data["ub"]
Q == Q' || warn("Q is not symmetric.")

costTrue = JLD.load("./MAT_FILES/objVals.jld")["nameValDict"][lowercase(filename)]

m = size(A,1)
n = size(Q,1)

Pa, qa, Aa, ba, K, typ = Converter.convertProblem(data)
println("Problem type: $(typ == 1 ? "Ax=b" : "l<=Ax<=u")")
settings = OSSDPSettings(rho=100.,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,scaling = 10,checkTermination = 1,eps_abs = 1e-3,eps_rel=1e-4,timelimit=10)
res,ws = OSSDP.solve(Pa,qa,Aa,ba,K,settings)

# # solve with OSQP
Aa2 = [A;eye(n)]
P2 = Q
q2 = c
if norm(rl-ru,Inf) < 1e-4
  u = [rl;ub]
  l = [rl;lb]
else
  l = [rl;lb]
  u = [ru;ub]
end

m2 = OSQP.Model()
OSQP.setup!(m2; P=P2, q=q2[:], A=Aa2, l=l[:], u=u[:],scaling=10,check_termination=1,verbose=false,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
res2 = OSQP.solve!(m2)


@testset "Maros Meszaros' QP Test Problems" begin

  @test abs(res.cost - costTrue) < 1e-2
  @test abs(res2.info.obj_val - costTrue) < 1e-2
  @test norm(res.x[1:n] - res2.x,Inf) < 1e-2


end