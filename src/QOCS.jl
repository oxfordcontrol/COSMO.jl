#__precompile__()
module QOCS

const DefaultFloat = Float64
const DefaultInt   = Int64

#export MathOptInterfaceQOCS
export  assemble!
#         reset!,
#         optimize!,
#         warmStart!,
#         Settings,
#         Model,
#         Result,
#         Constraint,
#         AbstractConvexSet

using SparseArrays,LinearAlgebra

include("./algebra.jl")
include("./projections.jl")
include("./types.jl")               # some types still need tidying
include("./settings.jl")            # unmodified - revisit
include("./constraint.jl")          # unmodified - revisit
include("./parameters.jl")          # unmodified - revisit
include("./residuals.jl")           # unmodified - revisit
include("./scaling.jl")             # set scaling / E scaling is broken
include("./kkt.jl")                 # unmodified - revisit.  Add lin solver type
include("./infeasibility.jl")       # unmodified - revisit.  Redundancy with composite set
#include("./printing.jl")            # unmodified - revisit (not used or uses SEDUMI K?)
include("./setup.jl")               # unmodified - revisit (very short - consolidate?)
include("./solver.jl")              # unmodified - revisit
include("./interface.jl")           # unmodified - revisit

# Is this module required?
# include("./Helper.jl")

end #end module
