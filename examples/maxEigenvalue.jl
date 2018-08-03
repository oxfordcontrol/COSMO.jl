# SDP example to find the maximum eigenvalue of a symmetric matrix A

# from: https://www.cs.cmu.edu/afs/cs.cmu.edu/academic/class/15859-f11/www/notes/lecture12.pdf
# Suppose A has eigenvalues λ1 ≥ λ2 . . . ≥ λn. Then the matrix tI − A has
# eigenvalues t−λ1, t−λ2, . . . , t−λn. Note that tI −A is psd exactly when all these eigenvalues
# are non-negative, and this happens for values t ≥ λ1
# This immediately gives us that
# λ1 = min{t | s.t. tI − A ≥ 0}

workspace()
include("../src/QOCS.jl")

using Base.Test
using QOCS

nn = 10
rng = MersenneTwister(7232)

@testset "Max Eigenvalue Problem with Random Matrices" begin

  for iii = 1:nn
    # generate symmetric test matrix A
    r = rand(rng,2:30)
    A = Symmetric(randn(rng,r,r))

    # solve the dual problem
    c = -vec(A)
    Aa = [vec(speye(r))';-speye(r^2)]
    b = [1.;zeros(r^2)]
    K = Cone(1,0,[],[r^2])
    P = spzeros(r^2,r^2)
    settings = QOCS.Settings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-4,eps_rel=1e-4)
    res,nothing = QOCS.solve(P,c[:],Aa,b[:],K,settings)
    println("$(iii)/$(nn) completed! Size of A: $(r), Number of Iterations $(res.iter).")

    # true solution
    λMaxTrue = maximum(eig(A)[1])
    @test abs(res.ν[1]-λMaxTrue) < 1e-2
   end
end

