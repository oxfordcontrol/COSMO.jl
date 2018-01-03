 # workspace()

# # include pure julia osqp solver
include("../osqp_julia.jl")

using OSQPSolver
using Base.Test
# # include python interface for native osqp solver
using PyCall
 @pyimport osqp as osqp
 @pyimport scipy.sparse as spar
 @pyimport numpy as np
 @pyimport scipy.sparse as pysparse

jlmat2pymat(S::SparseMatrixCSC) = pysparse.csc_matrix((S.nzval, S.rowval .- 1, S.colptr .- 1), shape=size(S))
pymat2jlmat(S::PyObject) = SparseMatrixCSC(S[:m], S[:n], S[:indptr] .+ 1, S[:indices] .+ 1, S[:data])


# Define solver settings for both solvers
α = 1.6;
ρ = 1.0;
σ = 1e-6;
eps_prim_inf = 1e-5;
eps_dual_inf = 1e-5;
eps_rel = 1e-5;
eps_abs = 1e-5;
max_iter = 2500;
verbose = 1;
scaling = 0;
polish = false;
warm_start = false;

# Simple problem
# m = 50
# n = 100
# A  = sparse(randn(m,n))
# l = -rand(m,1) * 2
# u = +rand(m,1) * 2
# P = Symmetric(sprand(n,n,0.1))
# P = P + n*speye(n)
# P = sparse(P)
# q = randn(n,1)

# Lasso problem (range 10-200)
n = 10
m = n*100
# data matrix
Ad = sprandn(m,n,0.5)
# true solution vector v with 50% zero elements and rest N(0,1/n)
v = randn(n)
v = v/sqrt(n)
r = rand(0:1,n)
v = [v[i]*r[i] for i=1:length(r)]
# right hand side vector b
b = Ad*v + randn(m)*0.25
λMax = norm(Ad'*b,Inf)
λ = 1/5 * λMax
# create QP matrizes
P = sparse([zeros(n,n+m+n);zeros(m,n) 2*eye(m) zeros(m,n); zeros(n,n+m+n)])
q = [zeros(m+n,1);λ*ones(n,1)]
A = [Ad -eye(m) zeros(m,n); eye(n) zeros(n,m) -eye(n); eye(n) zeros(n,m) eye(n)]
l = [b; -Inf*ones(n); zeros(n)]
u = [b; zeros(n); Inf*ones(n)]

# Primal infeasible problem
# n = 50
# m = 500
# Pt = sprandn(n, n, 0.6)
# P = Pt' * Pt
# q = randn(n, 1)
# A = sprandn(m, n, 0.8)
# u = 3 + randn(m, 1)
# l = -3 + randn(m, 1)

# # # Make random problem primal infeasible
# nhalf = Int64(floor(n/2))
# A[nhalf, :] = A[nhalf + 1, :]
# l[nhalf] = u[nhalf + 1] + 10 * rand()
# u[nhalf] = l[nhalf] + 0.5


# # Dual infeasible problem
# P = sparse(diagm([4; 0]))
# q = [0.0; 2]
# A = sparse([1.0 1; -1 1])
# l = [-Inf; -Inf]
# u = [2.0; 3]
# settings = qpSettings(rho=1.0,verbose=true)
# res = solveOSQP(P,q,A,l,u,settings)
# nothing

settings = qpSettings(rho=ρ,sigma=σ,alpha=α,eps_abs=eps_abs,eps_rel=eps_rel,eps_prim_inf=eps_prim_inf,eps_dual_inf=eps_dual_inf,max_iter=max_iter,verbose=Bool(verbose))
 res = solveOSQP(P,q,A,l,u,settings)
 nothing;

# convert julia sparse matrix type into python sparse matrix type
P_Py = jlmat2pymat(P)
A_Py = jlmat2pymat(A)

prob = osqp.OSQP()
prob[:setup](P_Py, q, A_Py, l, u, alpha=α,rho=ρ,sigma=σ,eps_prim_inf=eps_prim_inf,eps_dual_inf=eps_dual_inf,eps_rel=eps_rel,eps_abs=eps_abs,scaling=scaling,max_iter=max_iter,verbose=verbose,polish=polish,warm_start=warm_start)
res_py = prob[:solve]()

# # Compare solutions between Julia and Python
# # TODO: Improve check conditions
@testset "Julia Python Comparison" begin
    @test norm(res.x - res_py[:x]) < settings.eps_abs*10
    @test norm(res.y - res_py[:y]) < settings.eps_abs*10
    @test norm(res.cost - res_py[:info][:obj_val]) < settings.eps_abs*10
end