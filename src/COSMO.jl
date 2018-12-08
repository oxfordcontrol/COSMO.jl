#__precompile__()
module COSMO

using SparseArrays, LinearAlgebra, SuiteSparse

#export MathOptInterfaceCOSMO
export  assemble!, optimize!, warmStart!, reset!

const DefaultFloat = Float64
const DefaultInt   = Int64



include("./algebra.jl")
include("./projections.jl")
include("./types.jl")               # some types still need tidying
include("./settings.jl")            # TODO: unmodified - revisit
include("./constraint.jl")          # TODO: unmodified - revisit
include("./parameters.jl")          # TODO: unmodified - revisit
include("./residuals.jl")           # TODO: unmodified - revisit
include("./scaling.jl")             # TODO: set scaling / E scaling is broken
include("./kkt.jl")                 # TODO: unmodified - revisit.  Add lin solver type
include("./infeasibility.jl")       # TODO: stylistic fixes needed
include("./printing.jl")            # TODO: unmodified - revisit
include("./setup.jl")               # TODO: unmodified - revisit (short - consolidate?)
include("./solver.jl")              # TODO: unmodified - revisit
include("./interface.jl")           # TODO: unmodified - revisit

end #end module
