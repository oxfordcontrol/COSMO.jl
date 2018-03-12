# Test script to check if the matrix quilibration works
# the output matrix M should have similar lp norm along the rows
# the condition number κ(M)=σMax(M)/σMin(M) should decrease
# the transformation should hold MNew = S*M*S
workspace()

include("../Solver.jl")

using FactCheck, Scaling, OSSDP, OSSDPTypes


# create dummy settings object
settings = sdpSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=500,verbose=true,checkTermination=15,scaling=10)

# use ill-conditioned example data
#P = [1e5 0 0; 0 1e-5 0; 0 0 10]
P = zeros(4,4)
A = [1. 2 3 2; 4 100 3 0 ; 1 0 0 2]
m = size(A,1)
n = size(P,1)
q = [1e5;1;0;0]
b = [100;-1;1e4]
K = cone(0,0,[],[4])
p = problem(P,q,A,b,K)
sm = scaleMatrices()

# assemble matrix M
M = [P A';A zeros(m,m)]
# check condition number

κ = cond(full(M))
# check lp norm of rows of original matrix
rowNorms = [norm(M[i,:],2) for i in 1:size(M,1)]
colNorms = [norm(M[:,i],2) for i in 1:size(M,2)]
deltaRowNorm = maximum(rowNorms) - minimum(rowNorms)

# perform scaling

  D,E = scaleProblemSCS!(p,sm,settings)
S = [D zeros(n,m);zeros(m,n) E]
# reassamble scaled matrix
MNew = [p.P p.A';p.A zeros(m,m)]

# calculate condition number
κNew = cond(full(MNew))

# calculate lp norm of rows of scaled matrix
rowNormsNew = [norm(MNew[i,:],Inf) for i in 1:size(MNew,1)]
deltaRowNormNew = maximum(rowNormsNew) - minimum(rowNormsNew)


facts("Check matrix equilibration routine") do
  @fact maximum(abs.(MNew-(S*M*S)))--> less_than(1e-9)
  @fact κNew - κ --> less_than(0)
  @fact deltaRowNormNew - deltaRowNorm --> less_than(0)
  @fact isposdef(D) --> true
  @fact isposdef(E) --> true
end