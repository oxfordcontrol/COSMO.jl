#__precompile__()
module QOCS

const DefaultFloat = Float64
const DefaultInt   = Int64

#export MathOptInterfaceQOCS
# export  optimize!,
#         reset!,
#         assemble!,
#         warmStart!,
#         Settings,
#         Model,
#         Result,
#         Constraint,
#         AbstractConvexSet

using SparseArrays,LinearAlgebra

include("./algebra.jl")
include("./projections.jl")
include("./types.jl")
include("./constraint.jl")
include("./parameters.jl")
include("./residuals.jl")


include("./scaling.jl")

#old styles
# include("./Helper.jl")
# include("./KKT.jl")
#

# include("./Parameters.jl")
# include("./Infeasibility.jl")
# include("./Printing.jl")
# include("./Setup.jl")
# include("./Solver.jl")
# include("./Interface.jl")

end #end module
