# script to test new elimination tree from graph algorithm

workspace()
include("../../../src/Graph.jl")
include("../../../src/Tree.jl")

using GraphModule, TreeModule, Base.Test

A = [1 1 0 0 0; 1 1 1 0 0; 0 1 1 1 1; 0 0 1 1 1; 0 0 1 1 1]

g = Graph(sparse(A))


rng = Random.MersenneTwister(123554);
# Define number of test matrices
nn = 1000


function parentsFromOldAlgo(g)
  tic()
  t = TreeModule.createTreeFromGraph(g);
  tOld = toq()

  parOld = zeros(1:length(g.adjacencyList))
  for iii=1:length(t.nodes)
    n=t.nodes[iii]
    parOld[iii] = n.parent
  end
  return parOld, tOld
end


tNew3Arr = zeros(nn)
tNew4Arr = zeros(nn)
tOldArr = zeros(nn)

@testset "Elimination Tree Derivation from random Matrices" begin

    for iii=1:nn

        # take random dimension
        dim = rand(rng,50:200)
        density = rand(rng,0.1:0.1:0.6)
        # create random sparse matrix
        A = sprand(rng,dim,dim,density)
        A = A+A'

        # create graph from A and make Graph chordal
        g = Graph(A)
        # create Tree from Graph
        tic()
        par3 = TreeModule.etree3(g)
        tNew3 = toq()

        tic()
        par4 = TreeModule.etree4(g)
        tNew4 = toq()


        parOld,tOld = parentsFromOldAlgo(g)

        tNew3Arr[iii] = tNew3
        tNew4Arr[iii] = tNew4
        tOldArr[iii] = tOld


        @test par4 == parOld
        @test par3 == parOld
        println("$(iii) / $(nn)")
    end
end








