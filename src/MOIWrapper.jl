# This interface was derived from the OSQP MathOptInterface written by T. Koolen and B. Stellato
# https://github.com/oxfordcontrol/OSQP.jl/blob/master/src/MathOptInterfaceOSQP.jl


using MathOptInterface

export Optimizer

 MOI = MathOptInterface
 MOIU = MOI.Utilities
 CI = MOI.ConstraintIndex
 VI = MOI.VariableIndex

SparseTriplets = Tuple{Vector{Int}, Vector{Int}, Vector{<:Any}}

SingleVariable = MOI.SingleVariable
Affine = MOI.ScalarAffineFunction{Float64}
Quadratic = MOI.ScalarQuadraticFunction{Float64}
AffineConvertible = Union{Affine, SingleVariable}
VectorOfVariables = MOI.VectorOfVariables
VectorAffine = MOI.VectorAffineFunction{Float64}

Interval = MOI.Interval{Float64}
LessThan = MOI.LessThan{Float64}
GreaterThan = MOI.GreaterThan{Float64}
EqualTo = MOI.EqualTo{Float64}
IntervalConvertible = Union{LessThan, GreaterThan, EqualTo, Interval}


Zeros = MOI.Zeros
 #Nonnegatives = MOI.Nonnegatives
Nonpositives = MOI.Nonpositives
SOC = MOI.SecondOrderCone
PSD = Union{MOI.PositiveSemidefiniteConeSquare,MOI.PositiveSemidefiniteConeTriangle}
SupportedVectorSets = Union{Zeros, MOI.Nonnegatives, Nonpositives,SOC,PSD}

export printIdxmap,sortSets,assign_constraint_row_ranges!, processconstraints,constraint_rows, processobjective,processlinearterms!, symmetrize!, processconstraints!, constant, processconstant!, processlinearpart!, processconstraintset!

dimension(s::MOI.AbstractSet) = MOI.dimension(s)
dimension(::MOI.AbstractScalarSet) = 1
lower(::Zeros, i::Int) = 0.0
lower(::Nonnegatives, i::Int) = 0.0
lower(::Nonpositives, i::Int) = -Inf
upper(::Zeros, i::Int) = 0.0
upper(::Nonnegatives, i::Int) = Inf
upper(::Nonpositives, i::Int) = 0.0

# TODO: just use ∈ on 0.7 (allocates on 0.6):
function _contains(haystack, needle)
    for x in haystack
        x == needle && return true
    end
    false
end


##############################
# MAIN INTERFACE OBJECTS AND FUNCTIONS
##############################

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::COSMO.Model
    hasresults::Bool
    results::COSMO.Result{Float64}
    is_empty::Bool
    settings::COSMO.Settings
    sense::MOI.OptimizationSense
    objconstant::Float64
    constrconstant::Vector{Float64}
    rowranges::Dict{Int, UnitRange{Int}}

    function Optimizer()
        inner = COSMO.Model()
        hasresults = false
        results = COSMO.Result{Float64}()
        is_empty = true
        settings = COSMO.Settings()
        sense = MOI.MinSense
        objconstant = 0.
        constrconstant = Float64[]
        rowranges = Dict{Int, UnitRange{Int}}()
        new(inner, hasresults, results, is_empty, settings, sense, objconstant, constrconstant, rowranges)
    end
end

function Base.show(io::IO, obj::Optimizer)
    if obj.is_empty
        print(io,"Empty COSMO - Optimizer")
    elseif obj.hasresults
        print(io,"COSMO - Optimizer\n- Has results: $(obj.hasresults)\n- Objective constant: $(obj.objconstant)\n- Problem status: $(obj.results.status)\n- Optimal objective: $(round(obj.results.objVal, digits = 3))\n- Sense: $(obj.sense)\n- Iterations: $(obj.results.iter)\n- Solve time: $(round.(obj.results.times.solverTime*1000,digits=2))ms")
    else
        print(io,"COSMO - Optimizer\n- Has results: $(obj.hasresults)\n- Objective constant: $(obj.objconstant)\n- Sense: $(obj.sense)")
    end
end


hasresults(optimizer::Optimizer) = optimizer.hasresults

# MG: function to reset otimizer
function MOI.empty!(optimizer::Optimizer)
    optimizer.inner = COSMO.Model()
    optimizer.hasresults = false
    optimizer.results = COSMO.Result{Float64}()
    optimizer.is_empty = true
    optimizer.sense = MOI.MinSense # model parameter, so needs to be reset
    optimizer.objconstant = 0.
    optimizer.constrconstant = Float64[]
    empty!(optimizer.rowranges)
    optimizer
end

# MG: check if optimizer object is empty
MOI.is_empty(optimizer::Optimizer) = optimizer.is_empty


struct UnsupportedObjectiveError <: Exception end


##############################
# MODEL --> SOLVER FORMAT FUNCTIONS
##############################

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; copy_names=false)
    copy_names && error("Copying names is not supported.")
    MOI.empty!(dest)
    idxmap = MOIU.IndexMap(dest, src)
    assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
    A,b,convexSets = processconstraints(src, idxmap, dest.rowranges)
    COSMO.set!(dest.inner, P, q, A, b, convexSets)
    dest.is_empty = false
    idxmap
end

function MOI.optimize!(optimizer::Optimizer)
    optimizer.results = COSMO.optimize!(optimizer.inner, optimizer.settings)
    optimizer.hasresults = true
    nothing
end

# function add_constant(result::COSMO.Result{T}, constant::Real) where {T}
#     return new_result = COSMO.Result{T}(result.x, result.y, result.s, result.objVal + constant, result.iter, result.status, result.info, result.times);
# end

# create function for your solver that creates a index map, defined like this:
# struct IndexMap
#     varmap::Dict{MOI.VariableIndex, MOI.VariableIndex}
#     conmap::Dict{MOI.ConstraintIndex, MOI.ConstraintIndex}
# end
#FIX-ME: The order is important in my case, since s is ordered -> {0}, R+, SOCP, PSD
function MOIU.IndexMap(dest::Optimizer, src::MOI.ModelLike)
    idxmap = MOIU.IndexMap()
    # map model variable indices to solver variable indices
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = MOI.VariableIndex(i)
    end
    # map model constraint indices to solver constraint indices. For now this is important since solver expects following
    # order: {0}-variables, R+-variables, SOCP cones, psd cones
    LOCs = MOI.get(src, MOI.ListOfConstraints())
    #sort!(LOCs, by=x-> sortSets(x[2]))
    i = 0
    for (F, S) in LOCs
        MOI.supports_constraint(dest, F, S) || throw(UnsupportedConstraintError(F, S))
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = CI{F, S}(i)
        end
    end
    idxmap
end

# # returns a number for each supported set to ensure compatible order
# function sortSets(t::Type{<:MOI.AbstractSet})
#     (t <: MOI.EqualTo || t <: MOI.Zeros ) && return 1
#     (t <:MOI.GreaterThan || t <: MOI.LessThan || t <: MOI.Interval || t <: MOI.Nonnegatives || t <: MOI.Nonpositives) && return 2
#     (t <:MOI.SecondOrderCone) && return 3
#     (t <: MOI.PositiveSemidefiniteConeSquare) && return 4
#     return 5
# end

# Print idxmap for debugging purposes
function printIdxmap(idxmap::MOIU.IndexMap)
    println(">>Variable Map with $(length(idxmap.varmap)) entries:")
    dkeys = collect(keys(idxmap.varmap))
    dvalues = collect(values(idxmap.varmap))
    for i=1:length(dkeys)
        println("i=$(i): $(dkeys[i].value) => $(dvalues[i].value)")
    end

    println(">>Constraint Map with $(length(idxmap.conmap)) entries:")
    dkeys = collect(keys(idxmap.conmap))
    dvalues = collect(values(idxmap.conmap))
    for i=1:length(dkeys)
        println("i=$(i): $(dkeys[i].value) => $(dvalues[i].value)")
    end
end

# returns a Dictionary that maps a constraint index to a range
function assign_constraint_row_ranges!(rowranges::Dict{Int, UnitRange{Int}}, idxmap::MOIU.IndexMap, src::MOI.ModelLike)
    startrow = 1
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        # Returns Array of constraint indexes that match F,S, each constraint index is just a type-safe wrapper for Int64
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci_src in cis_src
            set = MOI.get(src, MOI.ConstraintSet(), ci_src)
            ci_dest = idxmap[ci_src]
            endrow = startrow + MOI.dimension(set) - 1
            rowranges[ci_dest.value] = startrow : endrow
            startrow = endrow + 1
        end
    end
end

function constraint_rows(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractScalarSet})
    rowrange = rowranges[ci.value]
    length(rowrange) == 1 || error()
    first(rowrange)
end
constraint_rows(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractVectorSet}) = rowranges[ci.value]
constraint_rows(optimizer::Optimizer, ci::CI) = constraint_rows(optimizer.rowranges, ci)


function processobjective(src::MOI.ModelLike, idxmap)
    sense = MOI.get(src, MOI.ObjectiveSense())
    n = MOI.get(src, MOI.NumberOfVariables())
    q = zeros(n)
    if sense != MOI.FeasibilitySense
        # FIXME: this is a hacky way of finding the objective function type and should be changed
        if typeof(src) <: MOIU.UniversalFallback
            obj_type = MOI.get(src.model,MOI.ObjectiveFunctionType())
        else
            obj_type = MOI.get(src, MOI.ObjectiveFunctionType())
        end

        if obj_type == MOI.SingleVariable
            fsingle = MOI.get(src, MOI.ObjectiveFunction{MOI.SingleVariable}())
            P = spzeros(n, n)
            q[idxmap[fsingle.variable].value] = 1
            c = 0.
        elseif obj_type == MOI.ScalarAffineFunction{Float64}
            faffine = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
            P = spzeros(n, n)
            # sets q
            processlinearterms!(q, faffine.terms, idxmap)
            c = faffine.constant
        elseif obj_type == MOI.ScalarQuadraticFunction{Float64}
            fquadratic = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}())
            I = [Int(idxmap[term.variable_index_1].value) for term in fquadratic.quadratic_terms]
            J = [Int(idxmap[term.variable_index_2].value) for term in fquadratic.quadratic_terms]
            V = [term.coefficient for term in fquadratic.quadratic_terms]
            symmetrize!(I, J, V)
            P = sparse(I, J, V, n, n)
            processlinearterms!(q, fquadratic.affine_terms, idxmap)
            c = fquadratic.constant
        else
            throw(UnsupportedObjectiveError())
        end
        sense == MOI.MaxSense && (rmul!(P, -1); rmul!(q, -1); c = -c)
    else
        P = spzeros(n, n)
        q = zeros(n)
        c = 0.
    end
    sense, P, q, c
end



function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm}, idxmapfun::Function = identity)
    for term in terms
        var = term.variable_index
        coeff = term.coefficient
        q[idxmapfun(var).value] += coeff
    end
end

function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm}, idxmap::MOIU.IndexMap)
    processlinearterms!(q, terms, var -> idxmap[var])
end


function symmetrize!(I::Vector{Int}, J::Vector{Int}, V::Vector)
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

function processconstraints(src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}})

    m = mapreduce(length, +, values(rowranges), init=0)

    b = zeros(Float64, m)
    constant = zeros(Float64, m)
    I = Int[]
    J = Int[]
    V = Float64[]
    COSMOconvexSets = Array{COSMO.AbstractConvexSet{Float64}}(undef, 0)
    # loop over constraints and modify A, l, u and constants
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        processconstraints!((I, J, V), b, COSMOconvexSets, constant, src, idxmap, rowranges, F, S)
    end
    # subtract constant from right hand side
    n = MOI.get(src, MOI.NumberOfVariables())
    A = sparse(I, J, V, m, n)
    return (A, b, COSMOconvexSets)
end

function processconstraints!(triplets::SparseTriplets, b::Vector, COSMOconvexSets::Array{COSMO.AbstractConvexSet{Float64}}, constant::Vector{Float64},
        src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}},
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    # loop over all constraints of same (F,S)
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        # handle cases where the sign of rows of A an b has to be flipped
        isa(s, Union{MOI.GreaterThan, MOI.Nonnegatives, SOC, PSD}) ? FLIP_SIGN = true : FLIP_SIGN = false

        rows = constraint_rows(rowranges, idxmap[ci])
        processConstant!(b, rows, f)
        processConstraint!(triplets, f, rows, idxmap, FLIP_SIGN)
        processSet!(b, rows, COSMOconvexSets, s)
    end
    nothing
end

##############################
# PROCESS FUNCTIONS
##############################
constant(f::MOI.SingleVariable) = 0
constant(f::MOI.ScalarAffineFunction) = f.constant


# process constant for functions Union{Affine, SingleVariable}
function processConstant!(b, row::Int, f::AffineConvertible)
    b[row] = constant(f)
    nothing
end

function processConstant!(b, rows::UnitRange{Int}, f::VectorOfVariables)
    b[rows] .= 0
    nothing
end

# process constant for functions VectorAffineFunction{Float64}
function processConstant!(b, rows::UnitRange{Int}, f::VectorAffine)
    for (i, row) in enumerate(rows)
        b[row] = f.constants[i]
    end
    nothing
end

# process function like f(x)= x1
function processConstraint!(triplets::SparseTriplets, f::MOI.SingleVariable, row::Int, idxmap, FLIP_SIGN::Bool)
    (I, J, V) = triplets
    col = idxmap[f.variable].value
    push!(I, row)
    push!(J, col)
    FLIP_SIGN ? push!(V, -1) : push!(V,1)
    nothing
end

# process function like f(x) = coeff*x_1
function processConstraint!(triplets::SparseTriplets, f::MOI.ScalarAffineFunction, row::Int, idxmap, FLIP_SIGN::Bool)
    (I, J, V) = triplets
    for term in f.terms
        var = term.variable_index
        coeff = term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        FLIP_SIGN ? push!(V, -coeff) : push!(V, coeff)
    end
end

function processConstraint!(triplets::SparseTriplets, f::MOI.VectorOfVariables, rows::UnitRange{Int}, idxmap, FLIP_SIGN::Bool)
    (I, J, V) = triplets
    for (i,var) in enumerate(f.variables)
        row = rows[i]
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        FLIP_SIGN ? push!(V, -1) : push!(V, 1)
    end
    nothing
end

# process function like f(x) = a'*x where x[1:n]
function processConstraint!(triplets::SparseTriplets, f::MOI.VectorAffineFunction, rows::UnitRange{Int}, idxmap, FLIP_SIGN::Bool)
    (I, J, V) = triplets
    for term in f.terms
        row = rows[term.output_index]
        var = term.scalar_term.variable_index
        coeff = term.scalar_term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        FLIP_SIGN ? push!(V, -coeff) : push!(V, coeff)
    end
end

# Special case for PositiveSemidefiniteConeTriangle, since this involves adding more rows to A and b
# function processConstraint!(triplets::SparseTriplets,f::MOI.VectorAffineFunction, rows::UnitRange{Int}, idxmap,FLIP_SIGN::Bool)
#     (I, J, V) = triplets
#     for term in f.terms
#         row = rows[term.output_index]
#         var = term.scalar_term.variable_index
#         coeff = term.scalar_term.coefficient
#         col = idxmap[var].value
#         push!(I, row)
#         push!(J, col)
#         FLIP_SIGN ? push!(V, -coeff) : push!(V, coeff)
#     end
# end

##############################
# PROCESS SETS
##############################

# process the following sets Union{LessThan, GreaterThan, EqualTo}
function processSet!(b::Vector, row::Int, cs, s::LessThan)
    b[row] += s.upper
    set = COSMO.Nonnegatives{Float64}(1)
    push!(cs, set)
    nothing
end
function processSet!(b::Vector, row::Int, cs, s::GreaterThan)
    b[row] -= s.lower
    set = COSMO.Nonnegatives{Float64}(1)
    push!(cs, set)
    nothing
end
function processSet!(b::Vector, row::Int, cs, s::EqualTo)
    b[row] += s.value
    set = COSMO.ZeroSet{Float64}(1)
    push!(cs, set)
    nothing
end

# process the following sets Zeros
function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::Zeros)
    set = COSMO.ZeroSet{Float64}(length(rows))
    push!(cs, set)
    nothing
end

# process the following sets Union{Zeros, Nonnegatives, Nonpositives}
function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.Nonnegatives)
    set = COSMO.Nonnegatives{Float64}(length(rows))
    push!(cs, set)
    nothing
end
function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::SOC)
    set = COSMO.SecondOrderCone{Float64}(length(rows))
    push!(cs, set)
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.PositiveSemidefiniteConeSquare)
    set = COSMO.PsdCone{Float64}(length(rows))
    push!(cs, set)
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.PositiveSemidefiniteConeTriangle)
    set = COSMO.PsdConeTriangle{Float64}(length(rows))
    push!(cs, set)
    nothing
end


##############################
# CAN GET / GET / SET : OPTIMIZER ATTRIBUTES
##############################

# MOI.canset(optimizer::Optimizer, a::COSMOAttribute) = isupdatable(a) || MOI.isempty(optimizer)
# function MOI.set!(optimizer::Optimizer, a::COSMOAttribute, value)
#     MOI.canset(optimizer, a) || error()
#     setting = Symbol(a)
#     setfield!(optimizer.settings,setting,value)
# end


## Supported Optimizer attributes, which  objective functions are supported:
# allow returning the inner model
MOI.get(optimizer::Optimizer, ::MOI.RawSolver) = optimizer.inner
# allow returning the number of results available
MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1


# Get objective function
function MOI.get(optimizer::Optimizer, a::MOI.ObjectiveValue)
    rawobj = optimizer.results.obj_val
    return ifelse(optimizer.sense == MOI.MaxSense, -rawobj, rawobj)
end

# Get solve time
MOI.get(optimizer::Optimizer, a::MOI.SolveTime) = optimizer.results.times.solverTime

MOI.get(optimizer::Optimizer, a::MOI.SolverName) = "COSMO"

# Get Termination Status
function MOI.get(optimizer::Optimizer, a::MOI.TerminationStatus)
    # Note that the :Dual_infeasible and :Primal_infeasible are mapped to MOI.Success
    # because COSMO can return a proof of infeasibility. For the same reason,
    # :Primal_infeasible_inaccurate is mapped to MOI.AlmostSuccess
    oscpstatus = optimizer.results.status
    if oscpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif oscpstatus == :Dual_infeasible
        MOI.Success
    elseif oscpstatus == :Primal_infeasible
        MOI.Success
    elseif oscpstatus == :Max_iter_reached
        MOI.IterationLimit
    elseif oscpstatus == :Solved
        MOI.Success
    end
end

# Get Primal Status
function MOI.get(optimizer::Optimizer, a::MOI.PrimalStatus)
    oscpstatus = optimizer.results.status
    if oscpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif oscpstatus == :Primal_infeasible
        MOI.InfeasibilityCertificate
    elseif oscpstatus == :Solved
        MOI.FeasiblePoint
    elseif oscpstatus == :Dual_infeasible
        MOI.InfeasibilityCertificate
    else
        MOI.UnknownResultStatus
    end
end

# Get Dual Status
function MOI.get(optimizer::Optimizer, a::MOI.DualStatus)
    oscpstatus = optimizer.results.info.status
    if oscpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif oscpstatus == :Dual_infeasible
        MOI.InfeasibilityCertificate
    elseif oscpstatus == :Primal_infeasible
        MOI.InfeasibilityCertificate
    elseif oscpstatus == :Solved
        MOI.FeasiblePoint
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate, :Non_convex (TODO: good idea? use COSMO.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end


## Variable attributes:


function MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::VI)
    return optimizer.results.x[vi.value]
end

function MOI.set(optimizer::Optimizer, a::MOI.VariablePrimalStart, vi::VI, value)
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    optimizer.warmstartcache.x[vi.value] = value
end


## Constraints:
function MOI.isvalid(optimizer::Optimizer, ci::CI)
    MOI.is_empty(optimizer) && return false
    ci.value ∈ keys(optimizer.rowranges)
end

function MOI.set(optimizer::Optimizer, a::MOI.ConstraintDualStart, ci::CI, value)
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    rows = constraint_rows(oscpstatusptimizer, ci)
    for (i, row) in enumerate(rows)
        optimizer.warmstartcache.y[row] = -value[i] # opposite dual convention
    end
    nothing
end


##############################
# SUPPORT OBJECTIVE AND SET FUNCTIONS
##############################



# allow objective functions that have the following template types
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{Quadratic}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true


## Supported Constraints
MOI.supports_constraint(optimizer::Optimizer, ::SingleVariable, ::GreaterThan) = true
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:AffineConvertible}, ::Type{<:IntervalConvertible}) = true
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:VectorAffine}, ::Type{<:SupportedVectorSets}) = true
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Union{VectorOfVariables,VectorAffine}}, ::Type{SOC}) = true
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Union{VectorOfVariables,VectorAffine}}, ::Type{<:PSD}) = true



