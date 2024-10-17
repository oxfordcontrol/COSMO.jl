# This interface was partially derived from the OSQP MathOptInterface written by T. Koolen and B. Stellato
# https://github.com/oxfordcontrol/OSQP.jl/blob/master/src/MathOptInterfaceOSQP.jl
# certain utility function are taken from the SCS  MathOptInterface:
# https://github.com/JuliaOpt/SCS.jl
import MathOptInterface as MOI
const MOIU = MOI.Utilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const SparseTriplets{Tv} = Tuple{Vector{<:Integer}, Vector{<:Integer}, Vector{Tv}}

const Affine = MOI.ScalarAffineFunction{<: AbstractFloat}
const Quadratic = MOI.ScalarQuadraticFunction{<: AbstractFloat}
const VectorAffine = MOI.VectorAffineFunction{<: AbstractFloat}

const Interval = MOI.Interval{<: AbstractFloat}
const GreaterThan = MOI.GreaterThan{<: AbstractFloat}
const IntervalConvertible = Union{GreaterThan, Interval}

const Zeros = MOI.Zeros
const SOC = MOI.SecondOrderCone
const PSDT = MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}
const ComplexPSDT = MOI.Scaled{MOI.HermitianPositiveSemidefiniteConeTriangle}
const SupportedVectorSets = Union{Zeros, MOI.Nonnegatives, SOC, PSDT, ComplexPSDT, MOI.ExponentialCone, MOI.DualExponentialCone, MOI.PowerCone, MOI.DualPowerCone}
const AggregatedSets = Union{Zeros, MOI.Nonnegatives, MOI.GreaterThan}
aggregate_equivalent(::Type{<: MOI.Zeros}) = COSMO.ZeroSet
aggregate_equivalent(::Type{<: Union{MOI.GreaterThan, MOI.Nonnegatives}}) = COSMO.Nonnegatives


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
    set_constant::Vector{T}
    constr_constant::Vector{T}
    rowranges::Dict{Int, UnitRange{Int}}
    idxmap::MOIU.IndexMap

    function Optimizer{T}(; user_settings...) where {T <: AbstractFloat}
        inner = COSMO.Model{T}()
        set_MOI_default_settings!(inner.settings)
        length(user_settings) > 0 && set_user_settings!(inner, user_settings)
        hasresults = false
        results = COSMO.Result{T}()
        is_empty = true
        sense = MOI.MIN_SENSE
        objconstant = zero(T)
        set_constant = T[]
        constr_constant = T[]
        rowranges = Dict{Int, UnitRange{Int}}()
        idxmap = MOIU.IndexMap()
        new(inner, hasresults, results, is_empty, sense, objconstant, set_constant, constr_constant, rowranges, idxmap)
    end
end
Optimizer(args...; kwargs...) = Optimizer{DefaultFloat}(args...; kwargs...)

function printIdxmap(idxmap::MOIU.IndexMap)
    println(">>Variable Map with $(length(idxmap.var_map)) entries:")
    dkeys = collect(keys(idxmap.var_map))
    dvalues = collect(values(idxmap.var_map))
    for i=1:length(dkeys)
        println("i=$(i): $(dkeys[i].value) => $(dvalues[i].value)")
    end

     println(">>Constraint Map with $(length(idxmap.con_map)) entries:")
    dkeys = collect(keys(idxmap.con_map))
    dvalues = collect(values(idxmap.con_map))
    for i=1:length(dkeys)
        println("i=$(i): $(dkeys[i].value) => $(dvalues[i].value)")
    end
end

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
    optimizer.set_constant = T[]
    optimizer.constr_constant = T[]   # the 5 in (3 * x1 + 2 * x2 + 5 ≤ 10 )
    optimizer.idxmap = MOIU.IndexMap()
    optimizer.rowranges = Dict{Int, UnitRange{Int}}()
    optimizer
end


# MG: check if optimizer object is empty
MOI.is_empty(optimizer::Optimizer) = optimizer.is_empty

##############################
# MODEL --> SOLVER FORMAT FUNCTIONS
##############################

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; copy_names = false)
    copy_names && error("Copying names is not supported.")
    MOI.empty!(dest)
    idxmap = MOIU.IndexMap(dest, src)
    LOCs = assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    processconstraints!(dest, src, idxmap, LOCs, dest.rowranges)
    pre_allocate_variables!(dest.inner)

    dest.is_empty = false
    dest.inner.states.IS_ASSEMBLED = true
    dest.idxmap = idxmap
    # pass attributes, e.g. objective, warm strarting vars, etc.
    dest.sense = MOI.get(src, MOI.ObjectiveSense())
    COSMO.pass_attributes!(dest, src, idxmap)

    # if no cost function is provided, allocate P, q with correct sizes
    dest.sense == MOI.FEASIBILITY_SENSE && COSMO.allocate_cost_variables!(dest)
    return idxmap
end
function MOI.optimize!(optimizer::Optimizer)
    optimizer.results = COSMO.optimize!(optimizer.inner)
    optimizer.hasresults = true
    nothing
end


function MOIU.IndexMap(dest::Optimizer, src::MOI.ModelLike)
    idxmap = MOIU.IndexMap()
    # map model variable indices to solver variable indices
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = MOI.VariableIndex(i)
    end
    # map model constraint indices to solver constraint indices. For now this is important since solver expects following
    # order: {0}-variables, R+-variables, SOCP cones, psd cones
    # LOCs = MOI.get(src, MOI.ListOfConstraintTypesPresent())
    # sort!(LOCs, by=x-> sort_sets(x[2]))
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        MOI.supports_constraint(dest, F, S) || throw(MOI.UnsupportedConstraint{F, S})
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = CI{F, S}(i)
        end
    end
    idxmap
end

# This is important for set merging
sort_sets(s::Type{<: MOI.Zeros}) = 1
sort_sets(s::Union{Type{ <: MOI.GreaterThan}, Type{ <: MOI.Nonnegatives}}) = 2
sort_sets(s::Type{MOI.Interval}) = 3
sort_sets(s::Type{<: MOI.AbstractSet}) = 4

# this has to be in the right order
# returns a Dictionary that maps a constraint index to a range
function assign_constraint_row_ranges!(rowranges::Dict{Int, UnitRange{Int}}, idxmap::MOIU.IndexMap, src::MOI.ModelLike)
    startrow = 1
    LOCs = MOI.get(src, MOI.ListOfConstraintTypesPresent())
    sort!(LOCs, by = x -> sort_sets(x[2]))
    for (F, S) in LOCs
        # Returns Array of constraint indices that match F,S, each constraint index is just a type-safe wrapper for Int
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci_src in cis_src
            set = MOI.get(src, MOI.ConstraintSet(), ci_src)
            ci_dest = idxmap[ci_src]
            endrow = startrow + MOI.dimension(set) - 1
            rowranges[ci_dest.value] = startrow : endrow
            startrow = endrow + 1
        end
    end
    return LOCs
end

function constraint_rows(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractScalarSet})
    rowrange = rowranges[ci.value]
    length(rowrange) == 1 || error()
    first(rowrange)
end
constraint_rows(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractVectorSet}) = rowranges[ci.value]

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
#     idxmap = dest.idxmap
#     n = MOI.get(dest, MOI.NumberOfVariables())
#     dest.inner.p.q = zeros(T, n)
#     dest.inner.p.q[idxmap[f.variable].value] = one(T)
#     dest.inner.p.P = spzeros(T, n, n)
#     dest.objconstant = zero(T)
#     dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
#     return nothing
# end

function process_objective!(dest::Optimizer{T}, f::MOI.ScalarAffineFunction{T}) where {T <: AbstractFloat}
    idxmap = dest.idxmap
    n = MOI.get(dest, MOI.NumberOfVariables())
    dest.inner.p.P = spzeros(T, n, n)

    dest.inner.p.q = zeros(T, n)
    processlinearterms!(dest.inner.p.q, f.terms, idxmap)

    dest.objconstant = f.constant
    dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
    return nothing
end

function process_objective!(dest::Optimizer{T}, f::MOI.ScalarQuadraticFunction{T}) where {T <: AbstractFloat}
    idxmap = dest.idxmap
    n = MOI.get(dest, MOI.NumberOfVariables())

    I = [Int(idxmap[term.variable_1].value) for term in f.quadratic_terms]
    J = [Int(idxmap[term.variable_2].value) for term in f.quadratic_terms]
    V = [term.coefficient for term in f.quadratic_terms]
    symmetrize!(I, J, V)
    dest.inner.p.P = sparse(I, J, V, n, n)

    dest.inner.p.q = zeros(T, n)
    processlinearterms!(dest.inner.p.q, f.affine_terms, idxmap)
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

function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm}, idxmap::MOIU.IndexMap)
    for term in terms
        var = term.variable
        coeff = term.coefficient
        q[idxmap[var].value] += coeff
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

function processconstraints!(optimizer::Optimizer{T}, src::MOI.ModelLike, idxmap, LOCs, rowranges::Dict{Int, UnitRange{Int}}) where {T <: AbstractFloat}
    
    if length(LOCs) > 0
        m = mapreduce(length, +, values(rowranges), init = 0)
    else
        # handle case of a problem without constraints
        m = 0
    end
    b = zeros(T, m)
    constant = zeros(T, m)
    I = Int[]
    J = Int[]
    V = T[]
    convex_sets = Array{COSMO.AbstractConvexSet{T}}(undef, 0)
    set_constant = zeros(T, m)
    
    if m > 0
        # loop over constraints and modify A, l, u and constants
        for (F, S) in LOCs
            processconstraints!((I, J, V), b, convex_sets, constant, set_constant, src, idxmap, rowranges, F, S)
        end
    else
        convex_sets = [COSMO.ZeroSet{T}(0)] #fall back if no constraints present
    end
    optimizer.set_constant = set_constant
    # subtract constant from right hand side
    n = MOI.get(src, MOI.NumberOfVariables())
    A = sparse(I, J, V, m, n)

    # store constraint matrices in inner model
    model = optimizer.inner
    model.p.A = A
    model.p.b = b
    model.p.C = CompositeConvexSet{T}(deepcopy(convex_sets))
    model.p.model_size = [size(A, 1); size(A,2 )]
    optimizer.constr_constant = constant
    return nothing
end

function processconstraints!(triplets::SparseTriplets{T}, b::Vector{T}, COSMOconvexSets::Array{COSMO.AbstractConvexSet{T}}, constant::Vector{T}, set_constant::Vector{T},
        src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}},
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet}) where {T <: AbstractFloat}
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    aggregate_dim = 0
    
    # loop over all constraints of same (F, S)
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        rows = constraint_rows(rowranges, idxmap[ci])
        aggregate_dim += length(rows)
        b[rows] = MOI.constant(f, T)
        constant[rows] = b[rows]
        if typeof(s) <: MOI.AbstractScalarSet && !(typeof(s) <: MOI.Interval)
            set_constant[rows] = MOI.constant(s)
        end
        processConstraint!(triplets, f, rows, idxmap, s)
        processSet!(b, rows, COSMOconvexSets, s)
    end

    process_aggregate_set!(b, COSMOconvexSets, aggregate_dim, cis_src, src, S)

    nothing
end

##############################
# PROCESS FUNCTIONS
##############################

# process function like f(x) = coeff*x_1
function processConstraint!(triplets::SparseTriplets{T}, f::MOI.ScalarAffineFunction, row::Int, idxmap::MOIU.IndexMap, s::MOI.AbstractSet) where {T <: AbstractFloat}
    (I, J, V) = triplets
    for term in f.terms
        push!(I, row)
        push!(J, idxmap[term.variable].value)
        push!(V, -term.coefficient)
    end
end

function processConstraint!(triplets::SparseTriplets{T}, f::MOI.VectorAffineFunction, rows::UnitRange{Int}, idxmap::MOIU.IndexMap, s::MOI.AbstractSet) where {T <: AbstractFloat}
    (I, J, V) = triplets

    vis_src = get_var_index.(f.terms)
    vis_dest = map(vi -> idxmap[vi], vis_src)

    A = sparse(output_index.(f.terms), map(x-> x.value, vis_dest), coefficient.(f.terms))
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    A_I, A_J, A_V = findnz(A)

    offset = first(rows) - 1

    append!(J, A_J)
    append!(V, -A_V)
    append!(I, A_I .+ offset)
end


##############################
# PROCESS SETS
##############################


process_aggregate_set!(b::Vector{T}, cs, dim::Int, cis_src, src, S::Type{<: MOI.AbstractSet}) where {T <: AbstractFloat} = false

function process_aggregate_set!(b::Vector{T}, cs, dim::Int, cis_src, src, S::Type{<: AggregatedSets}) where {T <: AbstractFloat}
   if length(cs) > 0 && cs[end] isa aggregate_equivalent(S){T}
        prev_dim = cs[end].dim    
        cs[end] = aggregate_equivalent(S){T}(dim + prev_dim)
   else
        push!(cs, aggregate_equivalent(S){T}(dim))
   end
end

function process_aggregate_set!(b::Vector{T}, cs, dim::Int, cis_src, src, s::Type{MOI.Interval{T}}) where {T <: AbstractFloat}
    l = zeros(T, length(cis_src))
    u = zeros(T, length(cis_src))
   for (k, ci) in enumerate(cis_src)
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        l[k] = s.lower
        u[k] = s.upper
    end
    push!(cs, COSMO.Box{T}(l, u))
    nothing
end


function processSet!(b::Vector{T}, row::Int, cs, s::GreaterThan) where {T <: AbstractFloat}
    b[row] -= s.lower
    # push!(cs, COSMO.Nonnegatives{T}(1))
    nothing
end

# do not process the AggregatedSets and Intervall constraints individually
function processSet!(b::Vector{T}, rows::Union{Int, UnitRange{Int}}, cs, s::Union{AggregatedSets, MOI.Interval{T}}) where {T <: AbstractFloat}
    nothing
end

function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::SOC) where {T <: AbstractFloat}
    push!(cs, COSMO.SecondOrderCone{T}(length(rows)))
    nothing
end


function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::PSDT) where {T <: AbstractFloat}
    push!(cs, COSMO.PsdConeTriangle{T, T}(length(rows)))
    nothing
end

function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::ComplexPSDT) where {T <: AbstractFloat}
    push!(cs, COSMO.PsdConeTriangle{T, Complex{T}}(length(rows)))
    nothing
end

function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::MOI.ExponentialCone) where {T <: AbstractFloat}
    push!(cs, COSMO.ExponentialCone{T}())
    nothing
end

function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::MOI.DualExponentialCone) where {T <: AbstractFloat}
    push!(cs, COSMO.DualExponentialCone{T}())
    nothing
end

function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::MOI.PowerCone) where {T <: AbstractFloat}
    push!(cs, COSMO.PowerCone{T}(T(s.exponent)))
    nothing
end

function processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::MOI.DualPowerCone) where {T <: AbstractFloat}
    push!(cs, COSMO.DualPowerCone{T}(T(s.exponent)))
    nothing
end




function pass_attributes!(dest::Optimizer{T}, src::MOI.ModelLike, idxmap::MOIU.IndexMap) where {T <: AbstractFloat}
    

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
            process_warm_start!(dest, MOI.VariablePrimalStart(), idxmap[vi], value)
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
                    process_warm_start!(dest, attr, idxmap[ci], value)
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
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Affine}, ::Type{<: IntervalConvertible}) = true
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:VectorAffine}, ::Type{<:SupportedVectorSets}) = true
