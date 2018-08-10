# Test chordal extension
workspace()
include("qdldl.jl")
include("/Users/Micha/Dropbox/Research/OSSDP/Code/src/Graph.jl")

using QDLDL,GraphModule, Base.Test



rng = Random.MersenneTwister(131123)



numProblems = 10





# g = Graph(A,false)
# println("Original Graph is chordal?: $(GraphModule.isPerfectOrdering(g))")
# # A = A  +(2*5+1)*speye(5);
# F = QDLDL.qdldl(A,logical=true,perm=nothing)
# L = F.L

# C = full(L+L') + eye(5)

# g_c = Graph(sparse(C))
# println("Extended Graph is chordal?: $(GraphModule.isPerfectOrdering(g_c))")

@testset "QDLDL routine for chordal embedding" begin

  for iii=1:numProblems

    dim = rand(rng,50:1:200);
    density=rand(rng,0.1:0.1:0.4)
    A = sprand(rng,dim,dim,density)

    A = 0.5*(A + A')
    A = A + (2*dim+1)*speye(dim);

    g = Graph(A,false)

    if (GraphModule.isPerfectOrdering(g) == true)
      continue
    else
      F = QDLDL.qdldl(A,logical=true,perm=nothing)
      L = F.L
      C = L+L' + speye(dim)
      g_c = Graph(C)
    end

    # check results
    @test GraphModule.isPerfectOrdering(g_c)

  end
end
