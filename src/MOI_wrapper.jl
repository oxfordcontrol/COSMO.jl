# This interface was partially derived from the OSQP MathOptInterface written by T. Koolen and B. Stellato
# https://github.com/oxfordcontrol/OSQP.jl/blob/master/src/MathOptInterfaceOSQP.jl
# certain utility function are taken from the SCS  MathOptInterface:
# https://github.com/JuliaOpt/SCS.jl
import MathOptInterface as MOI
const MOIU = MOI.Utilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const Affine = MOI.ScalarAffineFunction{<: AbstractFloat}
const Quadratic = MOI.ScalarQuadraticFunction{<: AbstractFloat}
const VectorAffine = MOI.VectorAffineFunction{<: AbstractFloat}

const SUPPORTED_CONES{T} = Union{
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.HyperRectangle{T},
    MOI.SecondOrderCone,
    MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
    MOI.Scaled{MOI.HermitianPositiveSemidefiniteConeTriangle},
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
    MOI.PowerCone{T},
    MOI.DualPowerCone{T},
}

# Mimicking SCS._SetConstants
struct _SetConstants{T}
    b::Vector{T}
    box_lower_or_power::Vector{T}
    box_upper::Vector{T}
    _SetConstants{T}() where {T} = new{T}(T[], T[], T[])
end

function Base.empty!(x::_SetConstants)
    empty!(x.b)
    empty!(x.box_lower_or_power)
    empty!(x.box_upper)
    return x
end

function Base.resize!(x::_SetConstants, n)
    resize!(x.b, n)
    resize!(x.box_lower_or_power, n)
    resize!(x.box_upper, n)
end

function MOI.Utilities.load_constants(x::_SetConstants, offset, f)
    MOI.Utilities.load_constants(x.b, offset, f)
    return
end

function MOI.Utilities.load_constants(
    x::_SetConstants{T},
    offset,
    set::MOI.HyperRectangle{T},
) where {T}
    I = offset .+ (eachindex(set.lower))
    x.box_lower_or_power[I] = set.lower
    x.box_upper[I] = set.upper
    return
end

function MOI.Utilities.load_constants(
    x::_SetConstants{T},
    offset,
    set::Union{MOI.PowerCone{T},MOI.DualPowerCone{T}},
) where {T}
    x.box_lower_or_power[offset+1] = set.exponent
    return
end

function MOI.Utilities.function_constants(x::_SetConstants, rows)
    return MOI.Utilities.function_constants(x.b, rows)
end

function MOI.Utilities.set_from_constants(x::_SetConstants, S, rows)
    return MOI.Utilities.set_from_constants(x.b, S, rows)
end

function MOI.Utilities.set_from_constants(
    x::_SetConstants{T},
    ::Type{S},
    rows,
) where {T,S<:Union{MOI.PowerCone{T},MOI.DualPowerCone{T}}}
    @assert length(rows) == 3
    return S(x.box_lower_or_power[first(rows)])
end

function MOI.Utilities.set_from_constants(
    x::_SetConstants{T},
    ::Type{MOI.HyperRectangle{T}},
    rows,
) where {T}
    return MOI.HyperRectangle(
        x.box_lower_or_power[rows],
        x.box_upper[rows],
    )
end

MOI.Utilities.@product_of_sets(
    Cones,
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.HyperRectangle{T},
    MOI.SecondOrderCone,
    MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
    MOI.Scaled{MOI.HermitianPositiveSemidefiniteConeTriangle},
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
    MOI.PowerCone{T},
    MOI.DualPowerCone{T},
)

const OptimizerCache{T} = MOI.Utilities.GenericModel{
    T,
    MOI.Utilities.ObjectiveContainer{T},
    MOI.Utilities.VariablesContainer{T},
    MOI.Utilities.MatrixOfConstraints{
        T,
        MOI.Utilities.MutableSparseMatrixCSC{
            T,
            Int,
            MOI.Utilities.OneBasedIndexing,
        },
        _SetConstants{T},
        Cones{T},
    },
}

##############################
# MAIN INTERFACE OBJECTS AND FUNCTIONS
##############################

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    inner::COSMO.Model{T}
    hasresults::Bool
    results::COSMO.Result{T}
    is_empty::Bool
    sense::MOI.OptimizationSense
    objconstant::T
    rowranges::Dict{Int, UnitRange{Int}}

    function Optimizer{T}(; user_settings...) where {T <: AbstractFloat}
        inner = COSMO.Model{T}()
        set_MOI_default_settings!(inner.settings)
        length(user_settings) > 0 && set_user_settings!(inner, user_settings)
        hasresults = false
        results = COSMO.Result{T}()
        is_empty = true
        sense = MOI.MIN_SENSE
        objconstant = zero(T)
        new(inner, hasresults, results, is_empty, sense, objconstant)
    end
end
Optimizer(args...; kwargs...) = Optimizer{DefaultFloat}(args...; kwargs...)

function Base.show(io::IO, obj::Optimizer)
    if obj.is_empty
        print(io,"Empty COSMO - Optimizer")
    elseif obj.hasresults
        print(io,"COSMO - Optimizer\n- Has results: $(obj.hasresults)\n- Objective constant: $(obj.objconstant)\n- Problem status: $(obj.results.status)\n- Optimal objective: $(round(obj.results.obj_val, digits = 3))\n- Sense: $(obj.sense)\n- Iterations: $(obj.results.iter)\n- Solve time: $(round.(obj.results.times.solver_time*1000,digits=2))ms")
    else
        print(io,"COSMO - Optimizer\n- Has results: $(obj.hasresults)\n- Objective constant: $(obj.objconstant)\n- Sense: $(obj.sense)")
    end
end

function set_user_settings!(inner::COSMO.Model, user_settings)
    for (k, v) in user_settings
        !in(k, fieldnames(typeof(inner.settings))) && error("The user setting $(k) is not defined.")
        setfield!(inner.settings, k, v)
    end
end

# Function to set MOI specific default values from COSMO that are different from plain COSMO
# the default behaviour for MOI solver is verbose printing turned on
function set_MOI_default_settings!(settings)
    settings.verbose = true
end


# MG: function to reset otimizer
function MOI.empty!(optimizer::Optimizer{T}) where {T <: AbstractFloat}
    COSMO.empty_model!(optimizer.inner)
    optimizer.hasresults = false
    optimizer.results = COSMO.Result{T}()
    optimizer.is_empty = true
    optimizer.sense = MOI.MIN_SENSE # model parameter, so needs to be reset
    optimizer.objconstant = zero(T)
    optimizer
end


# MG: check if optimizer object is empty
MOI.is_empty(optimizer::Optimizer) = optimizer.is_empty

##############################
# MODEL --> SOLVER FORMAT FUNCTIONS
##############################

function MOI.default_cache(::Optimizer{T}, ::Type{T}) where {T}
    return MOI.Utilities.UniversalFallback(OptimizerCache{T}())
end

function MOI.copy_to(
    dest::Optimizer{T},
    src::MOI.Utilities.UniversalFallback{OptimizerCache{T}},
) where {T}
    MOI.empty!(dest)

    processconstraints!(dest, src)
    pre_allocate_variables!(dest.inner)

    dest.is_empty = false
    dest.inner.states.IS_ASSEMBLED = true
    # pass attributes, e.g. objective, warm strarting vars, etc.
    dest.sense = MOI.get(src, MOI.ObjectiveSense())
    COSMO.pass_attributes!(dest, src)

    # if no cost function is provided, allocate P, q with correct sizes
    dest.sense == MOI.FEASIBILITY_SENSE && COSMO.allocate_cost_variables!(dest)
    return MOI.Utilities.identity_index_map(src)
end

function MOI.copy_to(dest::Optimizer{T}, src::MOI.ModelLike) where {T}
    cache = MOI.Utilities.UniversalFallback(OptimizerCache{T}())
    index_map = MOI.copy_to(cache, src)
    MOI.copy_to(dest, cache)
    return index_map
end

function MOI.optimize!(optimizer::Optimizer)
    optimizer.results = COSMO.optimize!(optimizer.inner)
    optimizer.hasresults = true
    nothing
end

# this is handled explicitly in `copy_to` before the objective function is processed
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
# MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveSense, sense) = nothing

function allocate_cost_variables!(optimizer::Optimizer)
    m, n = size(optimizer.inner.p.A)
    optimizer.inner.p.P = spzeros(n, n)
    optimizer.inner.p.q = zeros(n)
    return nothing
end

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{<:Union{Affine, Quadratic}}) = true
# function MOIU.load(dest::Optimizer{T}, ::MOI.ObjectiveFunction, f::MOI.VariableIndex) where {T <: AbstractFloat}
#     n = MOI.get(dest, MOI.NumberOfVariables())
#     dest.inner.p.q = zeros(T, n)
#     dest.inner.p.q[f.variable.value] = one(T)
#     dest.inner.p.P = spzeros(T, n, n)
#     dest.objconstant = zero(T)
#     dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
#     return nothing
# end

function process_objective!(dest::Optimizer{T}, f::MOI.ScalarAffineFunction{T}) where {T <: AbstractFloat}
    n = MOI.get(dest, MOI.NumberOfVariables())
    dest.inner.p.P = spzeros(T, n, n)

    dest.inner.p.q = zeros(T, n)
    processlinearterms!(dest.inner.p.q, f.terms)

    dest.objconstant = f.constant
    dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
    return nothing
end

function process_objective!(dest::Optimizer{T}, f::MOI.ScalarQuadraticFunction{T}) where {T <: AbstractFloat}
    n = MOI.get(dest, MOI.NumberOfVariables())

    I = [Int(term.variable_1.value) for term in f.quadratic_terms]
    J = [Int(term.variable_2.value) for term in f.quadratic_terms]
    V = [term.coefficient for term in f.quadratic_terms]
    symmetrize!(I, J, V)
    dest.inner.p.P = sparse(I, J, V, n, n)

    dest.inner.p.q = zeros(T, n)
    processlinearterms!(dest.inner.p.q, f.affine_terms)
    dest.objconstant = f.constant

    dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
    return nothing
end

function process_objective!(dest::Optimizer, f)
    throw(MOI.UnsupportedAttribute(f))
end

function apply_sense!(sense::MOI.OptimizationSense, P::AbstractMatrix{T}, q::AbstractVector{T}, c::T) where {T <: AbstractFloat}
    if sense == MOI.MAX_SENSE
        LinearAlgebra.rmul!(P, -one(T))
        LinearAlgebra.rmul!(q, -one(T))
        c = -c
    end
    return c
end

function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm})
    for term in terms
        var = term.variable
        coeff = term.coefficient
        q[var.value] += coeff
    end
end



function symmetrize!(I::Vector{Ti}, J::Vector{Ti}, V::Vector{Tv}) where {Tv <: AbstractFloat, Ti <: Integer}
    n = length(V)
    (length(I) == length(J) == n) || error()
    for i = 1 : n
        if I[i] != J[i]
            push!(I, J[i])
            push!(J, I[i])
            push!(V, V[i])
        end
    end
end

# The following utility functions are from https://github.com/JuliaOpt/SCS.jl

output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
get_var_index(t::MOI.ScalarAffineTerm) = t.variable
get_var_index(t::MOI.VectorAffineTerm) = get_var_index(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)

function processconstraints!(optimizer::Optimizer{T}, src::MOI.ModelLike) where {T <: AbstractFloat}
    convex_sets = COSMO.AbstractConvexSet{T}[]

    for (F, S) in MOI.get(cone_Ab, MOI.ListOfConstraintTypesPresent())
        process_sets!(convex_sets, cone_Ab, F, S)
    end
    
    if isempty(convex_sets)
        push!(convex_sets, COSMO.ZeroSet{T}(0)) #fall back if no constraints present
    end

    # store constraint matrices in inner model
    model = optimizer.inner
    model.p.A = -convert(SparseMatrixCSC{T,Int}, src.constraints.coefficients)
    model.p.b = src.constraints.constants.b
    model.p.C = CompositeConvexSet{T}(deepcopy(convex_sets))
    model.p.model_size = [size(A, 1), size(A, 2)]
    return nothing
end

function processSets!(convex_sets, src, ::Type{MOI.VectorAffineFunction{T}}, ::Type{MOI.Zeros}) where {T}
    push!(convex_sets, COSMO.ZeroSet{T}(MOI.Utilities.num_rows(src, MOI.Zeros)))
    return nothing
end

function processSets!(convex_sets, src, ::Type{MOI.VectorAffineFunction{T}}, ::Type{MOI.Nonnegatives}) where {T}
    push!(convex_sets, COSMO.Nonnegatives{T}(MOI.Utilities.num_rows(src, MOI.Nonnegatives)))
    return nothing
end

function processSets!(convex_sets, src, ::Type{MOI.VectorAffineFunction{T}}, ::Type{MOI.HyperRectangle{T}}) where {T}
    offset = MOIU.num_rows(src, MOI.Zeros) + MOIU.num_rows(src, MOI.Nonnegatives)
    rows = offset .+ (1:MOIU.num_rows(src, MOI.HyperRectangle{T}))
    push!(convex_sets, COSMO.Box{T}(src.constants.box_lower_or_power[rows], src.constants.box_upper[rows]))
    return nothing
end

function processSets!(convex_sets, src, ::Type{MOI.VectorAffineFunction{T}}, ::Type{S}) where {T, S}
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        set = MOI.get(src, MOI.ConstraintSet(), ci)::S
        processSet!(convex_sets, set, T)
    end
    return nothing
end

function processSet!(cs, s::MOI.SecondOrderCone, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.SecondOrderCone{T}(MOI.dimension(s)))
    nothing
end

function processSet!(cs, s::MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.PsdConeTriangle{T, T}(MOI.dimension(s)))
    nothing
end

function processSet!(cs, s::MOI.Scaled{MOI.HermitianPositiveSemidefiniteConeTriangle}, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.PsdConeTriangle{T, Complex{T}}(MOI.dimension(s)))
    nothing
end

function processSet!(cs, ::MOI.ExponentialCone, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.ExponentialCone{T}())
    nothing
end

function processSet!(cs, ::MOI.DualExponentialCone, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.DualExponentialCone{T}())
    nothing
end

function processSet!(cs, s::MOI.PowerCone, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.PowerCone{T}(T(s.exponent)))
    nothing
end

function processSet!(cs, s::MOI.DualPowerCone, ::Type{T}) where {T <: AbstractFloat}
    push!(cs, COSMO.DualPowerCone{T}(T(s.exponent)))
    nothing
end



function pass_attributes!(dest::Optimizer{T}, src::MOI.ModelLike) where {T <: AbstractFloat}
    

    model_attributes = MOI.get(src, MOI.ListOfModelAttributesSet())
    for attr in model_attributes
        if attr != MOI.ObjectiveSense() && attr != MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}() && attr != MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}() && attr != MOI.Name() 
            throw(MOI.UnsupportedAttribute(attr))
        end
    end

    # copy objective function
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    obj = MOI.get(src, MOI.ObjectiveFunction{F}())
    process_objective!(dest, obj) 
    
    # copy variable attributes, e.g. VariablePrimalStart
    has_primal_start = false
    var_attr = MOI.get(src, MOI.ListOfVariableAttributesSet())
    # We don't support MOI.VariableNames at this point
    for attr in var_attr
        if attr isa MOI.VariablePrimalStart
            has_primal_start = true
        end
    end
    if has_primal_start
        vis_src = MOI.get(src, MOI.ListOfVariableIndices())
        for vi in vis_src 
            value = MOI.get(src, MOI.VariablePrimalStart(), vi)
            process_warm_start!(dest, MOI.VariablePrimalStart(), vi, value)
        end
    end

    # Copy constraint attributes, e.g. ConstraintPrimalStart, ConstraintDualStart
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F, S}())
            if attr == MOI.ConstraintName()
                # skip
            elseif attr == MOI.ConstraintPrimalStart() || attr == MOI.ConstraintDualStart()
                for ci in cis_src
                    value = MOI.get(src, attr, ci)
                    process_warm_start!(dest, attr, ci, value)
                end
            else
                throw(MOI.UnsupportedAttribute(attr))
            end
        end
    end
    return nothing
end


##############################
# GET / SET : OPTIMIZER ATTRIBUTES
##############################


## Supported Optimizer attributes, which  objective functions are supported:
# allow returning the inner model
MOI.get(optimizer::Optimizer, ::MOI.RawSolver) = optimizer.inner
# allow returning the number of results available
MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = optimizer.hasresults ? 1 : 0
MOI.get(optimizer::Optimizer, ::MOI.NumberOfVariables) = length(optimizer.inner.vars.x)
function MOI.get(optimizer::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i in 1:length(optimizer.inner.vars.x)]
end

# Get objective function
function MOI.get(optimizer::Optimizer, a::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, a)
    rawobj = optimizer.results.obj_val
    if optimizer.sense == MOI.MAX_SENSE
        return -rawobj - optimizer.objconstant
    else
        return rawobj + optimizer.objconstant
    end
end


MOI.get(optimizer::Optimizer, ::MOI.SolverName) = "COSMO"
MOI.get(optimizer::Optimizer, ::MOI.SolverVersion) = "v" * string(COSMO.version())
MOI.get(optimizer::Optimizer, ::MOI.NumberOfThreads) = Threads.nthreads()

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
   if optimizer.hasresults 
        return string(optimizer.results.status)
   else
        return ""
   end
end

function MOI.get(optimizer::Optimizer, a::MOI.SolveTimeSec) 
    if optimizer.hasresults
        return optimizer.results.times.solver_time
    else
        return NaN
    end
end


# Get Termination Status
function MOI.get(optimizer::Optimizer, a::MOI.TerminationStatus)
    # Note that the :Dual_infeasible and :Primal_infeasible are mapped to MOI.Success
    # because COSMO can return a proof of infeasibility. For the same reason,
    # :Primal_infeasible_inaccurate is mapped to MOI.AlmostSuccess
    optimizer.hasresults || return MOI.OPTIMIZE_NOT_CALLED

    status = optimizer.results.status

    if status == :Dual_infeasible
        return MOI.DUAL_INFEASIBLE
    elseif status == :Primal_infeasible
        return MOI.INFEASIBLE
    elseif status == :Max_iter_reached
        return MOI.ITERATION_LIMIT
    elseif status == :Solved
        return MOI.OPTIMAL
    elseif status == :Time_limit_reached
        return MOI.TIME_LIMIT
    else
        return MOI.NUMERICAL_ERROR
    end
end

# Get Primal Status
function MOI.get(optimizer::Optimizer, a::MOI.PrimalStatus)
    if !(1 <= a.result_index <= MOI.get(optimizer, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

    status = optimizer.results.status
    if status == :Primal_infeasible
        return MOI.INFEASIBLE_POINT
    elseif status == :Solved
        return MOI.FEASIBLE_POINT
    elseif status == :Max_iter_reached || status == :Time_limit_reached
        settings = optimizer.inner.settings
        if is_primal_feasible(optimizer.results.info, settings)
            return MOI.FEASIBLE_POINT
        elseif is_primal_nearly_feasible(optimizer.results.info, settings)
            return MOI.NEARLY_FEASIBLE_POINT
        else
            return MOI.INFEASIBLE_POINT
        end
    elseif status == :Dual_infeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        return MOI.NO_SOLUTION
    end
end

# Get Dual Status
function MOI.get(optimizer::Optimizer, a::MOI.DualStatus)
    if !(1 <= a.result_index <= MOI.get(optimizer, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

    status = optimizer.results.status
    if status == :Dual_infeasible
        return MOI.INFEASIBLE_POINT
    elseif status == :Primal_infeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif status == :Solved
        return MOI.FEASIBLE_POINT
    elseif status == :Max_iter_reached || status == :Time_limit_reached
        settings = optimizer.inner.settings
        if is_dual_feasible(optimizer.results.info, settings)
            return MOI.FEASIBLE_POINT
        elseif is_dual_nearly_feasible(optimizer.results.info, settings)
            return MOI.NEARLY_FEASIBLE_POINT
        else
            return MOI.INFEASIBLE_POINT
        end
    else
        return MOI.NO_SOLUTION
    end
end

struct RawResult <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::RawResult) = true
MOI.get(optimizer::Optimizer, ::RawResult) = optimizer.results
MOIU.map_indices(::Function, r::COSMO.Result) = r


_unshift(optimizer::Optimizer, offset, value, s) = value
_unshift(optimizer::Optimizer, offset, value, s::Type{<:MOI.AbstractScalarSet}) = value + optimizer.set_constant[offset]
_shift(optimizer::Optimizer, offset, value, s) = value
_shift(optimizer::Optimizer, offset, value, s::Type{<:MOI.AbstractScalarSet}) = value - optimizer.set_constant[offset]

function MOI.get(optimizer::Optimizer, a::MOI.ConstraintPrimal, ci_src::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, a)
    # offset = constroffset(optimizer, ci)
    # rows = constrrows(optimizer, ci)
    rows = COSMO.constraint_rows(optimizer.rowranges, ci_src)
    c_primal = optimizer.results.s[rows]
    # (typeof(c_primal) <: AbstractArray && length(c_primal) == 1) && (c_primal = first(c_primal))
    return _unshift(optimizer, rows, c_primal, S)
end

function MOI.get(optimizer::Optimizer, a::MOI.ConstraintDual, ci_src::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, a)
    rows = constraint_rows(optimizer.rowranges, ci_src)
    return optimizer.results.y[rows]
end

## Variable attributes:
function MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::VI)
    MOI.check_result_index_bounds(optimizer, a)
    return optimizer.results.x[vi.value]
end

## Warm starting:
MOI.supports(::Optimizer, a::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}) = true
function process_warm_start!(optimizer::Optimizer, a::MOI.VariablePrimalStart, vi::VI, value::Real)
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    COSMO.warm_start_primal!(optimizer.inner, value, vi.value)
end
process_warm_start!(optimizer::Optimizer, a::MOI.VariablePrimalStart, vi::VI, value::Nothing) = nothing

MOI.supports(::Optimizer, a::MOI.ConstraintPrimalStart, ::Type{<:MOI.ConstraintIndex}) = true
function process_warm_start!(optimizer::Optimizer, a::MOI.ConstraintPrimalStart, ci::CI{<:MOI.AbstractFunction, S}, value) where S <: MOI.AbstractSet
    (value == nothing || isa(value, Array{Nothing, 1})) && return nothing
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    rows = constraint_rows(optimizer.rowranges, ci)
    # this undoes the shifting that was used in get(MOI.ConstraintPrimal)
    value = _shift(optimizer, rows, value, S)
    COSMO.warm_start_slack!(optimizer.inner, value, rows)
end

MOI.supports(::Optimizer, a::MOI.ConstraintDualStart, ::Type{<:MOI.ConstraintIndex}) = true
function process_warm_start!(optimizer::Optimizer, a::MOI.ConstraintDualStart, ci::CI{<:MOI.AbstractFunction, S}, value) where S <: MOI.AbstractSet
    (value == nothing || isa(value, Array{Nothing, 1})) && return nothing
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    rows = constraint_rows(optimizer.rowranges, ci)
    COSMO.warm_start_dual!(optimizer.inner, value, rows)
    nothing
end


 MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true
function MOI.set(optimizer::Optimizer, p::MOI.RawOptimizerAttribute, value)
    setfield!(optimizer.inner.settings, Symbol(p.name), value)
end

function MOI.get(optimizer::Optimizer, p::MOI.RawOptimizerAttribute)
    if in(Symbol(p.name), fieldnames(typeof(optimizer.inner.settings)))
        return getfield(optimizer.inner.settings, Symbol(p.name))
    end
    error("RawOptimizerAttribute with name $(p.name) is not set.")
end


# MOI.Silent == COSMO.Settings.verbose = false
MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.inner.settings.verbose = !value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = !optimizer.inner.settings.verbose

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(optimizer::Optimizer{T}, ::MOI.TimeLimitSec, value::Real) where {T <: AbstractFloat}
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_limit"), T(value))
end

function MOI.set(optimizer::Optimizer{T}, attr::MOI.TimeLimitSec, ::Nothing) where {T <: AbstractFloat}
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_limit"), zero(T))
end
function MOI.get(optimizer::Optimizer, ::MOI.TimeLimitSec)
    return MOI.get(optimizer, MOI.RawOptimizerAttribute("time_limit"))
end

"""
    ADMMIterations()
The number of ADMM iterations completed during the solve.
"""
struct ADMMIterations <: MOI.AbstractModelAttribute
end

MOI.is_set_by_optimize(::ADMMIterations) = true

function MOI.get(optimizer::Optimizer, ::ADMMIterations)
    return optimizer.results.iter
end

##############################
# SUPPORT OBJECTIVE AND SET FUNCTIONS
##############################


## Supported Constraints
MOI.supports_constraint(optimizer::Optimizer{T}, ::Type{MOI.VectorAffineFunction{T}}, ::Type{<:SUPPORTED_CONES{T}}) where {T} = true
