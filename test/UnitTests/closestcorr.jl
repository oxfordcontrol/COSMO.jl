# Closest Correlation matrix test matrix


# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

using COSMO, Test, LinearAlgebra, SparseArrays, Random


# this function creates a matrix A that slices out the diagonal entries Xii of a vectorized square matrix x=vec(X)
function create_diagonal_extractor(n)
  A = spzeros(n, n^2)
  A[1,1 ] = 1
  for iii = 2:n-1
    col = (iii - 1) * (n + 1)
    A[iii, col + 1] = 1
  end
  A[n, n^2] = 1
  return A
end

function diagonal_equal_one(A, atol)
  for iii = 1:size(A, 1)
    (abs(A[iii, iii] - 1) > atol) && return false
  end
  return true
end

rng = Random.MersenneTwister(12345)
xMin = -1.
xMax = 1.
n = 50
C = xMin .+ randn(rng, n, n)*(xMax - xMin)
c = vec(C)

isposdef(C) && warn("The perturbed correlation matrix is still pos def.")


n2 = n^2
m = n + n2

P = sparse(1.0I, n2, n2)
q = -vec(C)
r = 0.5*vec(C)' * vec(C)

A1 = create_diagonal_extractor(n)
b1 = -ones(n)
cs1 = COSMO.Constraint(A1, b1, COSMO.ZeroSet)

A2 = sparse(1.0I, n2, n2)
b2 = zeros(n2)
cs2 = COSMO.Constraint(A2, b2, COSMO.PsdCone)
constraints = [cs1; cs2]


model = COSMO.Model()
assemble!(model, P, q, constraints)

res = COSMO.optimize!(model)

Xsol = reshape(res.x, n, n)

@testset "Closest Correlation Matrix - SDP Problems" begin
  @test res.status == :Solved
  @test diagonal_equal_one(Xsol, 1e-5)
  @test minimum(eigen(Xsol).values) > -1e-3

end
nothing
