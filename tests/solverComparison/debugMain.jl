workspace()
include("../../src/Solver.jl")
include("../meszaros/ConvertProblem.jl")
include("./Compare.jl")

using OSQP, OSSDP, Compare,MAT


problem = "QAFIRO"
dirPath = "../meszaros/bart_meszaros_data/"
gc()
# load problem data
data = matread("$(dirPath)"*"$(problem).mat")
P1, q1,r1, A1, b1, K1 = loadMeszarosData(data,"OSSDP")
P2, q2,r2, A2, l2, u2 = loadMeszarosData(data,"OSQP")
pDims = getMeszarosDim(data)