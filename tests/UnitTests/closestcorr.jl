# Closest Correlation matrix test matrix


# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

using OSSDP, Base.Test


# this function creates a matrix A that slices out the diagonal entries Xii of a vectorized square matrix x=vec(X)
function createDiagonalExtractor(n)
  A = spzeros(n,n^2)
  A[1,1] = 1
  for iii=2:n-1
    col = (iii-1)*(n+1)
    A[iii,col+1] = 1
  end
  A[n,n^2] = 1
  return A
end

function diagonalEqualOne(A,atol)
  for iii=1:size(A,1)
    (abs(A[iii,iii] - 1) > atol) && return false
  end
  return true
end

rng = MersenneTwister(12345)
xMin = -50.
xMax = 50.
n = 50
C = xMin+randn(rng,n,n)*(xMax-xMin)
c = vec(C)

isposdef(C) && warn("The perturbed correlation matrix is still pos def.")


n2 = n^2
m = n+n2

P = speye(n2)
q = -vec(C)
r = 0.5*vec(C)'*vec(C)
b = [ones(n);zeros(n2)]
A = createDiagonalExtractor(n)
Aa = [A; -speye(n2)]
# specify cone
Kf = n
Kl = 0
Kq = []
Ks = [n^2]

K = OSSDPTypes.Cone(Kf,Kl,Kq,Ks)
settings = OSSDPTypes.OSSDPSettings()

res,nothing = OSSDP.solve(P,q,Aa,b,K,settings);
Xsol = reshape(res.x,n,n)

@testset "Closest Correlation Matrix - SDP Problems" begin
  @test res.status == :solved
  @test diagonalEqualOne(Xsol,1e-5)
  @test minimum(eig(Xsol)[1]) > -1e-3

end


