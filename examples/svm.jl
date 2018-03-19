# Test script to test solver for a svm problem in qp format
workspace()
include("../src/Solver.jl")

using Base.Test
using OSSDP, OSSDPTypes, PyPlot, Helper

# SVM problem as QP
# min w'w + λ 1/n ∑ ϵ_i
# s.t. yi(wxi-b) ≥ 1- ϵ_i, ∀ i

rng = MersenneTwister(1234)

# create problem data
n = 2
m = 200
mHalf = Int(m/2)
y = [ones(mHalf);-1*ones(mHalf)]
X = zeros(m,n)
Xupp = (1/sqrt(n))*sprandn(rng,mHalf,n,0.5)
Xlow = (1/sqrt(n))*sprandn(rng,mHalf,n,0.5)
X[1:mHalf,:] = Xupp + duplicateSparsityPattern(Xupp)*(1/n)
X[mHalf+1:end,:] = Xlow - duplicateSparsityPattern(Xlow)*(1/n)


#  n = 2
#  m = 9
# # # y=[-1.;1]
# # # X = [1 2;1 0]
# # y=[-1.;1;1;-1]
# # X = [-1 2;2 2;2 -2;-1 0]
# # non separable example (might be stupid)
# y=[-1.;-1;-1;-1;-1;1;1;1;1]
# X = [-2 -2;-2 0;-2 2; -2 4; 15 0;2 -2;2 0; 2 2; 2 4]
λ = 1

# create augmented matrices for QP problem
A = [-diagm(y)*X y -eye(m);
     zeros(m,n+1) -eye(m)]
b = [ones(m);zeros(m)]
q = [zeros(n+1);1/m*ones(m)]
P = blkdiag(2*λ*speye(n),spzeros(1+m,1+m))

K = OSSDPTypes.Cone(0,2*m,[],[])
# define example problem
settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=true,checkTermination=1,scaling = 0)
 res,nothing = OSSDP.solve(P,q,A,b,K,settings)


# # solve with osqp
# A2 = sparse([A[:,1:n+1+m];zeros(m,n+1) eye(m)])
# q2 = q[1:n+1+m]
# P2 = sparse(P[1:n+1+m,1:n+1+m])
# l = [ones(m);zeros(m)]
# u = Inf*ones(2*m)
# m2 = OSQP.Model()
# OSQP.setup!(m2; P=P2, q=q2, A=A2, l=l, u=u,scaling=10,check_termination=1,verbose=true,adaptive_rho = false,eps_abs = 1e-3,eps_rel=1e-3)
# res2 = OSQP.solve!(m2)





w = res.x[1:n]
bi = res.x[n+1]
x1 = collect(Float64,-n:0.001:n)
x2 = (bi - (w[1].*x1))/w[2]
PyPlot.figure(1,facecolor="white",figsize=(12,5))
PyPlot.xticks(fontsize=15)
PyPlot.yticks(fontsize=15)
PyPlot.scatter(X[1:mHalf+1,1],X[1:mHalf+1,2],color="r")
PyPlot.hold(true)
PyPlot.scatter(X[mHalf+2:end,1],X[mHalf+2:end,2],color="b")
PyPlot.plot(x1,x2,color="k")
PyPlot.grid(true)
PyPlot.axis([-n n -n n]')

PyPlot.xlabel(L"x_1",fontsize=15)
PyPlot.ylabel(L"x_2",fontsize=15)
