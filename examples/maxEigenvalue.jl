# SDP example to find the maximum eigenvalue of a symmetric matrix A

# from: https://www.cs.cmu.edu/afs/cs.cmu.edu/academic/class/15859-f11/www/notes/lecture12.pdf
# Suppose A has eigenvalues λ1 ≥ λ2 . . . ≥ λn. Then the matrix tI − A has
# eigenvalues t−λ1, t−λ2, . . . , t−λn. Note that tI −A is psd exactly when all these eigenvalues
# are non-negative, and this happens for values t ≥ λ1
# This immediately gives us that
# λ1 = min{t | s.t. tI − A ≥ 0}

# we solve the dual problem
# max tr(A*X), s.t. tr(I*X) = 1, X ≥ 0,
# the max eigenvalue is then equal to the dual variable ν

workspace()
include("../src/Solver.jl")

using Base.Test
using OSSDP, OSSDPTypes

nn = 10
rng = MersenneTwister(7232)

@testset "Max Eigenvalue Problem with Random Matrices" begin

  for iii = 1:nn
    # generate symmetric test matrix A
    r = rand(rng,2:50)
    A = Symmetric(randn(rng,r,r))

    # solve the dual problem
    c = -vec(A)
    Aa = vec(eye(r))'
    b = [1.]
    K = Cone(0,0,[],[r^2])
    P = zeros(r^2,r^2)
    settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-4,eps_rel=1e-4)
    res,nothing = OSSDP.solve(P,c,Aa,b,K,settings)
    println("$(iii)/$(nn) completed! Size of A: $(r), Number of Iterations $(res.iter).")

    # true solution
    λMaxTrue = maximum(eig(A)[1])
    @test abs(res.ν[1]-λMaxTrue) < 1e-2
   end
end

