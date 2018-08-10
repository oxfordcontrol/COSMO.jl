# File to check and compare qdldl routine against built-in cholesky routine
include("qdldl.jl")

using QDLDL, Base.Test


nDense = 1000
nSparse = 1000
numProblems = nDense + nSparse

rng = Random.MersenneTwister(131123)


function generateSparsePosDefMatrix(rng,n::Int64,density::Float64)
      X = sprand(rng,n,n,density)
      X = 0.5(X+X')
      # make random matrix pos def
      X = X+(2*n+1)*speye(n);
      return X
  end


function generatePosDefMatrix(rng,n::Int64)
      X = rand(rng,n,n)
      Q, R = qr(X)
      eigs = rand(rng,n).*(5 .-0.1) .+ 0.1
      X = Q*Matrix(Diagonal(eigs))*Q'
      X = 0.5*(X+X')
      return X
  end


solveTimesQDLDL = zeros(numProblems)
solveTimesCHOL = zeros(numProblems)


@testset "QDLDL routine for random problems" begin

  for iii=1:numProblems

    dim = rand(rng,50:1:200);

    if iii <=nSparse
      density=rand(rng,0.1:0.1:0.4)
      A = generateSparsePosDefMatrix(rng,dim,density)
    else
      A = sparse(generatePosDefMatrix(rng,dim))
    end
    xtrue = rand(rng,dim)
    btrue = A*xtrue


    # solve with QDLDL
    tic()
    F = QDLDL.qdldl(A)
    x_qdldl = solve(F,btrue)
    solveTimesQDLDL[iii] = toq()

    # solve with cholesky
    tic()
    W = cholfact(A)
    x_chol = W\btrue
    solveTimesCHOL[iii] = toq()

    # check results
    @test norm(xtrue-x_qdldl,Inf) < 1e-9
    @test norm(x_chol-x_qdldl,Inf) < 1e-9

  end
end

# print mean solve times for both routines
println("\n\n$(nSparse) SPARSE PROBLEMS\n"*"-"^30*"\nMean solve time QLDL routine: \t$(mean(solveTimesQDLDL[1:nSparse])*1000)ms,\nMean solve Cholesky routine: \t$(mean(solveTimesCHOL[1:nSparse])*1000)ms.\n"*"-"^30)
println("\n$(nDense) DENSE PROBLEMS\n"*"-"^30*"\nMean solve time QLDL routine: \t$(mean(solveTimesQDLDL[nSparse+1:end])*1000)ms,\nMean solve Cholesky routine: \t$(mean(solveTimesCHOL[nSparse+1:end])*1000)ms.\n"*"-"^30)