# Closest Correlation matrix test matrix


# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

using COSMO, Test, LinearAlgebra, SparseArrays, Random

if @isdefined UnitTestFloat
    T = UnitTestFloat
else
    T = Float64
end
T = Float32


# this function creates a matrix A that slices out the diagonal entries Xii of a vectorized square matrix x=vec(X)
function create_diagonal_extractor(n, type::Type{T} = Float64) where {T <: AbstractFloat}
  A = spzeros(type, n, n^2)
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
xMin = -one(T)
xMax = one(T)
n = 50
C = xMin .+ randn(rng, T, n, n)*(xMax - xMin)
c = vec(C)


n2 = n^2
m = n + n2

P = spdiagm(0 => ones(T, n2))
q = -vec(C)
r = T(0.5) * dot(c, c)

A1 = create_diagonal_extractor(n, T)
b1 = -ones(T, n)
cs1 = COSMO.Constraint(A1, b1, COSMO.ZeroSet)

A2 = spdiagm(0 => ones(T, n2))
b2 = zeros(T, n2)
cs2 = COSMO.Constraint(A2, b2, COSMO.PsdCone)
constraints = [cs1; cs2]


model = COSMO.Model{T}()
assemble!(model, P, q, constraints)

res = COSMO.optimize!(model)

Xsol = reshape(res.x, n, n)

@testset "Closest Correlation Matrix - SDP Problems" begin
  @test res.status == :Solved
  @test diagonal_equal_one(Xsol, 1e-5)
  @test minimum(real.(eigen(Xsol).values)) > -1e-3

end
nothing
