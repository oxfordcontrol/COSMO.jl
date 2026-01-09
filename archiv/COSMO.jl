__precompile__()
module COSMO

using SparseArrays, LinearAlgebra, SuiteSparse, QDLDL, Pkg, DataStructures, Requires, Printf, IterTools
using Reexport
using COSMOAccelerators

using Reexport
@reexport using COSMOAccelerators

export assemble!, warm_start!, empty_model!, update!

const RealOrComplex{T <: Real} = Union{T, Complex{T}}
const DefaultFloat = Float64
const DefaultInt   = LinearAlgebra.BlasInt


include("./linear_solver/kktsolver.jl")
# optional dependencies
function __init__()
    @require Pardiso="46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2" include("./linear_solver/kktsolver_pardiso.jl")

    @require IterativeSolvers="42fd0dbc-a981-5370-80f2-aaf504508153" begin
        @require LinearMaps="7a12625a-238d-50fd-b39a-03d52299707e" include("./linear_solver/kktsolver_indirect.jl")
    end
end

function version()
    v"0.8.9"
end


include("./algebra.jl")
include("./projections.jl")
include("./chordal_decomposition/trees.jl")
include("./chordal_decomposition/clique_graph.jl")
include("./chordal_decomposition/clique_merging.jl")
include("./settings.jl")
include("./types.jl")
include("./accelerator_interface.jl")

include("./constraint.jl")
include("./parameters.jl")
include("./residuals.jl")
include("./scaling.jl")
include("./infeasibility.jl")
include("./chordal_decomposition/transformations.jl")
include("./chordal_decomposition/chordal_decomposition.jl")
include("./printing.jl")
include("./setup.jl")
include("./solver.jl")
include("./interface.jl")
include("./MOI_wrapper.jl")


end #end module
