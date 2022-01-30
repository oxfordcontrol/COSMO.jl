__precompile__()
module COSMO

using SparseArrays, LinearAlgebra, SuiteSparse, QDLDL, Pkg, DataStructures, Requires, Printf, IterTools, GenericLinearAlgebra, MutableArithmetics
using Reexport
using COSMOAccelerators

using Reexport
@reexport using COSMOAccelerators

export assemble!, warm_start!, empty_model!, update!

const DefaultFloat = Float64
const DefaultInt   = LinearAlgebra.BlasInt
const MA = MutableArithmetics


include("./kktsolver.jl")
# optional dependencies
function __init__()
    @require Pardiso="46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2" include("kktsolver_pardiso.jl")

    @require IterativeSolvers="42fd0dbc-a981-5370-80f2-aaf504508153" begin
        @require LinearMaps="7a12625a-238d-50fd-b39a-03d52299707e" include("./kktsolver_indirect.jl")
    end
end

function version()
    v"0.8.1"
end


include("./algebra.jl")
include("./projections.jl")
#include("./acceleration.jl")
include("./trees.jl")
include("./clique_graph.jl")
include("./clique_merging.jl")
include("./settings.jl")
include("./types.jl")
include("./accelerator_interface.jl")

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
include("./MOI_wrapper.jl")

# include("precompile.jl")
# _precompile_()

end #end module
