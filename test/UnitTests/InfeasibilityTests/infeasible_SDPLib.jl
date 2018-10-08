# Test of the infeasibility detection feature of the COSMO solver on problems from the SDPLibrary
# SDP problems are given in the following format
# problem: max  tr( F0*Y)
#          s.t. tr(Fi Y) = ci, for i=1,...,m
#               Y ⪴ 0

# RESULTS SO FAR:
# infp1 primal form --> detected as prim inf, dual form --> detected as dual inf
# infp2 primal form --> detected as prim inf, dual form --> detected as dual inf
# infd1 primal form --> not detected, dual form --> not detected (both forms dual inf according to MOSEK)
# infd2 primal form --> not detected, dual form --> not detected (both forms dual inf according to MOSEK)

# workspace()
# include("../../../src/Solver.jl")
# using JLD
# using Base.Test
# using COSMO
#  using JuMP, Mosek
using JLD

problems = ["infp1.jld";"infp2.jld";"infd1.jld";"infd2.jld"]

@testset "Infeasible problems from SDP Library" begin
  for kkk = 1:length(problems)
    data = JLD.load("../sdplib/"*problems[kkk])
    F = data["F"]
    c = data["c"]
    m = data["m"]
    n = data["n"]
    optVal = data["optVal"]
    kkk <= 2 ? (problem_type = :Primal_infeasible) : (problem_type=:Dual_infeasible)

    # Rewrite problem to COSMO compatible format:
    # min   1/2 x'Px + q'x
    # s.t.  Ax + s = b
    #       s in K
    # -----------------------
    # primal form
    # -----------------------
    P = zeros(m,m)
    q = vec(c)
    A = zeros(n^2,m)
    for iii = 1:m
      A[:,iii] = -vec(F[iii+1])
    end
    b = -vec(F[1])
    Kf = 0
    Kl = 0
    Kq = []
    Ks = [n^2]
    # model = Model(solver=MosekSolver())
    # @variable(model, x[1:m])
    # @variable(model, S[1:n,1:n],SDP)
    # @objective(model, Min, q'*x)
    # s = vec(S)
    # @constraint(model, A*x+s .== b)
    # status = JuMP.solve(model)

    # -----------------------
    # dual form
    # -----------------------
    # P = zeros(n^2,n^2)
    # q = -vec(F[1])
    # A = zeros(m,n^2)
    # for iii = 1:m
    #   A[iii,:] = vec(F[iii+1])'
    # end
    # A = [A; -speye(n^2)]
    # b = [c;zeros(n^2)]
    # Kf = m
    # Kl = 0
    # Kq = []
    # Ks = [n^2]

    # model = Model(solver=MosekSolver())
    # @variable(model, x[1:n^2])
    # @variable(model, s1[1:m])
    # @variable(model, S2[1:n,1:n],SDP)
    # s = [s1;vec(S2)]
    # @objective(model, Min, q'*x)
    # @constraint(model, A*x+s .== b)
    # status = JuMP.solve(model)
    K = Cone(Kf,Kl,Kq,Ks)
    settings = Settings(max_iter=3000,verbose=false,check_termination=1,checkInfeasibility=50,scaling = 10 ,scaleFunc=2,adaptive_rho=true,eps_abs=1e-4,eps_rel=1e-4)
    res = COSMO.solve(P,q,A,b,K,settings);
    @test res.status == problem_type
  end
end
