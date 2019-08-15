__precompile__()
module COSMO

using SparseArrays, LinearAlgebra, SuiteSparse, QDLDL, Pkg


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

include("./updatable_thin_qr.jl")
include("./algebra.jl")
include("./projections.jl")
include("./settings.jl")            # TODO: unmodified - revisit
include("./trees.jl")
include("./types.jl")               # some types still need tidying
include("./constraint.jl")          # TODO: unmodified - revisit
include("./parameters.jl")          # TODO: unmodified - revisit
include("./residuals.jl")           # TODO: unmodified - revisit
include("./scaling.jl")             # TODO: set scaling / E scaling is broken
include("./infeasibility.jl")       # TODO: stylistic fixes needed
include("./chordal_decomposition.jl")
include("./printing.jl")            # TODO: unmodified - revisit
include("./setup.jl")               # TODO: unmodified - revisit (short - consolidate?)
include("./solver.jl")              # TODO: unmodified - revisit
include("./interface.jl")           # TODO: unmodified - revisit
include("./MOIWrapper.jl")

export extract_upper_triangle, populate_upper_triangle, UpdatableQ, add_column!, add_columns!

end #end module
