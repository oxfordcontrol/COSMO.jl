function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(COSMO, Symbol("#19#20")) && precompile(Tuple{getfield(COSMO, Symbol("#19#20")),Tuple{SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.ExponentialCone{Float64}}})
    isdefined(COSMO, Symbol("#19#20")) && precompile(Tuple{getfield(COSMO, Symbol("#19#20")),Tuple{SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.Nonnegatives{Float64}}})
    isdefined(COSMO, Symbol("#19#20")) && precompile(Tuple{getfield(COSMO, Symbol("#19#20")),Tuple{SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.PsdConeTriangle{Float64}}})
    isdefined(COSMO, Symbol("#19#20")) && precompile(Tuple{getfield(COSMO, Symbol("#19#20")),Tuple{SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.SecondOrderCone{Float64}}})
    isdefined(COSMO, Symbol("#19#20")) && precompile(Tuple{getfield(COSMO, Symbol("#19#20")),Tuple{SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.ZeroSet{Float64}}})
    isdefined(COSMO, Symbol("#25#26")) && precompile(Tuple{getfield(COSMO, Symbol("#25#26")),Tuple{SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.SecondOrderCone{Float64}}})
    precompile(Tuple{Type{CholmodKKTSolver},SparseArrays.SparseMatrixCSC{Float64,Int64},SparseArrays.SparseMatrixCSC{Float64,Int64},Float64,Array{Float64,1}})
    precompile(Tuple{typeof(Base.Broadcast._broadcast_getindex_evalf),typeof(MathOptInterface.add_constraint),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.ExponentialCone})
    precompile(Tuple{typeof(Base.Broadcast._broadcast_getindex_evalf),typeof(MathOptInterface.add_constraint),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Nonnegatives})
    precompile(Tuple{typeof(Base.Broadcast._broadcast_getindex_evalf),typeof(MathOptInterface.add_constraint),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.PositiveSemidefiniteConeTriangle})
    precompile(Tuple{typeof(Base.Broadcast._broadcast_getindex_evalf),typeof(MathOptInterface.add_constraint),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.SecondOrderCone})
    precompile(Tuple{typeof(Base.Broadcast._broadcast_getindex_evalf),typeof(MathOptInterface.add_constraint),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Zeros})
    precompile(Tuple{typeof(Base.Broadcast.restart_copyto_nonleaf!),Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},S} where S,1},Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.SecondOrderCone},1},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1},Tuple{Base.OneTo{Int64}},typeof(MathOptInterface.add_constraint),Tuple{Base.RefValue{MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}}},Base.Broadcast.Extruded{Array{Union{MathOptInterface.SingleVariable, MathOptInterface.VectorAffineFunction{Float64}},1},Tuple{Bool},Tuple{Int64}},Base.Broadcast.Extruded{Array{Any,1},Tuple{Bool},Tuple{Int64}}}},MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Nonnegatives},Int64,Base.OneTo{Int64},Int64,Int64})
    precompile(Tuple{typeof(COSMO.admm_step!),Array{Float64,1},COSMO.SplitVector{Float64},Array{Float64,1},SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},Array{Float64,1},Array{Float64,1},Array{Float64,1},CholmodKKTSolver{Float64,Int64},Array{Float64,1},Array{Float64,1},Array{Float64,1},Float64,Float64,Int64,Int64,COSMO.CompositeConvexSet{Float64}})
    precompile(Tuple{typeof(COSMO.apply_sense!),MathOptInterface.OptimizationSense,SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1},Float64})
    precompile(Tuple{typeof(COSMO.dual_kkt_condition!),Array{Float64,1},Array{Float64,1},SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1},Array{Float64,1},SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1}})
    precompile(Tuple{typeof(COSMO.kkt_col_norms!),SparseArrays.SparseMatrixCSC{Float64,Int64},SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1},Array{Float64,1}})
    precompile(Tuple{typeof(COSMO.optimize!),COSMO.Workspace{Float64}})
    precompile(Tuple{typeof(COSMO.pass_attributes!),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Function})
    precompile(Tuple{typeof(COSMO.pre_allocate_variables!),COSMO.Workspace{Float64}})
    precompile(Tuple{typeof(COSMO.primal_kkt_condition!),Array{Float64,1},SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1},COSMO.SplitVector{Float64},Array{Float64,1}})
    precompile(Tuple{typeof(COSMO.processconstraints!),Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}},Array{Float64,1},Array{COSMO.AbstractConvexSet{Float64},1},Array{Float64,1},Array{Float64,1},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Dict{Int64,UnitRange{Int64}},Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.ExponentialCone}})
    precompile(Tuple{typeof(COSMO.processconstraints!),Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}},Array{Float64,1},Array{COSMO.AbstractConvexSet{Float64},1},Array{Float64,1},Array{Float64,1},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Dict{Int64,UnitRange{Int64}},Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.Nonnegatives}})
    precompile(Tuple{typeof(COSMO.processconstraints!),Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}},Array{Float64,1},Array{COSMO.AbstractConvexSet{Float64},1},Array{Float64,1},Array{Float64,1},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Dict{Int64,UnitRange{Int64}},Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.PositiveSemidefiniteConeTriangle}})
    precompile(Tuple{typeof(COSMO.processconstraints!),Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}},Array{Float64,1},Array{COSMO.AbstractConvexSet{Float64},1},Array{Float64,1},Array{Float64,1},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Dict{Int64,UnitRange{Int64}},Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.SecondOrderCone}})
    precompile(Tuple{typeof(COSMO.processconstraints!),Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}},Array{Float64,1},Array{COSMO.AbstractConvexSet{Float64},1},Array{Float64,1},Array{Float64,1},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Dict{Int64,UnitRange{Int64}},Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.Zeros}})
    precompile(Tuple{typeof(COSMO.rectify_scaling!),SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},COSMO.SecondOrderCone{Float64}})
    precompile(Tuple{typeof(COSMO.scale_data!),SparseArrays.SparseMatrixCSC{Float64,Int64},SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1},Array{Float64,1},Diagonal{Float64,Array{Float64,1}},Diagonal{Float64,Array{Float64,1}},Float64})
    precompile(Tuple{typeof(COSMO.scale_data!),SparseArrays.SparseMatrixCSC{Float64,Int64},SparseArrays.SparseMatrixCSC{Float64,Int64},Array{Float64,1},Array{Float64,1},LinearAlgebra.UniformScaling{Bool},Diagonal{Float64,Array{Float64,1}},Float64})
    precompile(Tuple{typeof(MathOptInterface.Bridges.Objective.add_all_bridges),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},Type})
    precompile(Tuple{typeof(MathOptInterface.Bridges.Variable.add_all_bridges),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},Type})
    precompile(Tuple{typeof(MathOptInterface.Utilities._pass_attributes),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},Bool,MathOptInterface.Utilities.IndexMap,Array{MathOptInterface.AbstractConstraintAttribute,1},Tuple{DataType},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.ExponentialCone},1}},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.ExponentialCone},1}},typeof(MathOptInterface.Utilities.load)})
    precompile(Tuple{typeof(MathOptInterface.Utilities._pass_attributes),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},Bool,MathOptInterface.Utilities.IndexMap,Array{MathOptInterface.AbstractConstraintAttribute,1},Tuple{DataType},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Nonnegatives},1}},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Nonnegatives},1}},typeof(MathOptInterface.Utilities.load)})
    precompile(Tuple{typeof(MathOptInterface.Utilities._pass_attributes),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},Bool,MathOptInterface.Utilities.IndexMap,Array{MathOptInterface.AbstractConstraintAttribute,1},Tuple{DataType},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.PositiveSemidefiniteConeTriangle},1}},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.PositiveSemidefiniteConeTriangle},1}},typeof(MathOptInterface.Utilities.load)})
    precompile(Tuple{typeof(MathOptInterface.Utilities._pass_attributes),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},Bool,MathOptInterface.Utilities.IndexMap,Array{MathOptInterface.AbstractConstraintAttribute,1},Tuple{DataType},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.SecondOrderCone},1}},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.SecondOrderCone},1}},typeof(MathOptInterface.Utilities.load)})
    precompile(Tuple{typeof(MathOptInterface.Utilities._pass_attributes),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},Bool,MathOptInterface.Utilities.IndexMap,Array{MathOptInterface.AbstractConstraintAttribute,1},Tuple{DataType},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Zeros},1}},Tuple{Array{MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Zeros},1}},typeof(MathOptInterface.Utilities.load)})
    precompile(Tuple{typeof(MathOptInterface.Utilities.load),Optimizer,MathOptInterface.ObjectiveFunction{MathOptInterface.SingleVariable},MathOptInterface.SingleVariable})
    precompile(Tuple{typeof(MathOptInterface.add_variables),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},Int64})
    precompile(Tuple{typeof(MathOptInterface.get),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.ConstraintDual,MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Nonnegatives}})
    precompile(Tuple{typeof(MathOptInterface.get),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.ConstraintDual,MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Zeros}})
    precompile(Tuple{typeof(MathOptInterface.get),MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}},MathOptInterface.ObjectiveValue})
    precompile(Tuple{typeof(copy),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1},Tuple{Base.OneTo{Int64}},typeof(MathOptInterface.add_constraint),Tuple{Base.RefValue{MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}}},Array{Union{MathOptInterface.SingleVariable, MathOptInterface.VectorAffineFunction{Float64}},1},Array{Any,1}}}})
    precompile(Tuple{typeof(empty_model!),COSMO.Workspace{Float64}})
    precompile(Tuple{typeof(setproperty!),COSMO.ProblemData{Float64},Symbol,SparseArrays.SparseMatrixCSC{Float64,Int64}})
end
