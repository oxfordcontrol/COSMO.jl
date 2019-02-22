# SDP example to find the maximum eigenvalue of a symmetric matrix A

# from: https://www.cs.cmu.edu/afs/cs.cmu.edu/academic/class/15859-f11/www/notes/lecture12.pdf
# Suppose A has eigenvalues λ1 ≥ λ2 . . . ≥ λn. Then the matrix tI − A has
# eigenvalues t−λ1, t−λ2, . . . , t−λn. Note that tI −A is psd exactly when all these eigenvalues
# are non-negative, and this happens for values t ≥ λ1
# This immediately gives us that
# λ1 = min{t | s.t. tI − A ≥ 0}
using Test
using COSMO, SparseArrays, LinearAlgebra, Random

nn = 10
rng = MersenneTwister(7232)

@testset "Max Eigenvalue Problem with Random Matrices" begin

  for iii = 1:nn
    # generate symmetric test matrix A
    r = rand(rng, 2:30)
    A = Symmetric(randn(rng, r, r))

    # solve the dual problem
    c = -vec(A)
    A1 = -vec(sparse(1.0I, r, r))'
    A2 = sparse(1.0I, r^2, r^2)
    b1 = 1.
    b2 = zeros(r^2)

    constraint1 = COSMO.Constraint(A1, b1, COSMO.ZeroSet)
    constraint2 = COSMO.Constraint(A2, b2, COSMO.PsdCone)
    P = spzeros(r^2, r^2)

    settings = COSMO.Settings(check_termination=1, scaling = 0)

    model = COSMO.Model()
    assemble!(model, P, c,[constraint1; constraint2], settings = settings)
    res = COSMO.optimize!(model);



    println("$(iii)/$(nn) completed! Size of A: $(r), Number of Iterations $(res.iter).")

    # true solution
    λMaxTrue = maximum(eigen(A).values)
    @test abs(res.y[1] - λMaxTrue) < 1e-2
   end
end
nothing
