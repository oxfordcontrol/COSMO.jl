__precompile__()
module COSMO

using SparseArrays, LinearAlgebra, SuiteSparse, QDLDL, Pkg, DataStructures


export assemble!, warmStart!, empty_model!

const DefaultFloat = Float64
const DefaultInt   = LinearAlgebra.BlasInt


include("./kktsolver.jl")
# optional dependencies
if in("Pardiso",keys(Pkg.installed()))
    include("./kktsolver_pardiso.jl")
end
if in("IterativeSolvers", keys(Pkg.installed())) && in("LinearMaps", keys(Pkg.installed()))
    include("./kktsolver_indirect.jl")
end

include("./algebra.jl")
include("./projections.jl")
include("./trees.jl")
include("./clique_merging.jl")
include("./settings.jl")
include("./types.jl")
include("./constraint.jl")
include("./parameters.jl")
include("./residuals.jl")
include("./scaling.jl")
include("./infeasibility.jl")
include("./transformations.jl")
include("./chordal_decomposition.jl")
include("./printing.jl")
include("./setup.jl")
include("./solver.jl")
include("./interface.jl")           
include("./MOIWrapper.jl")


end #end module
