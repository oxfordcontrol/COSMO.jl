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
PSDS = MOI.PositiveSemidefiniteConeSquare
PSDT = MOI.PositiveSemidefiniteConeTriangle
PSD = Union{MOI.PositiveSemidefiniteConeSquare,MOI.PositiveSemidefiniteConeTriangle}
SupportedVectorSets = Union{Zeros, MOI.Nonnegatives, Nonpositives,SOC, PSD}

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
    sense::MOI.OptimizationSense
    objconstant::Float64
    set_constant::Vector{Float64}
    constr_constant::Vector{Float64}
    rowranges::Dict{Int, UnitRange{Int}}
    idxmap::MOIU.IndexMap

    function Optimizer(; user_settings...)
        inner = COSMO.Model()
        length(user_settings) > 0 && set_user_settings!(inner, user_settings)
        hasresults = false
        results = COSMO.Result{Float64}()
        is_empty = true
        sense = MOI.MinSense
        objconstant = 0.
        set_constant = Float64[]
        constr_constant = Float64[]
        rowranges = Dict{Int, UnitRange{Int}}()
        idxmap = MOIU.IndexMap()
        new(inner, hasresults, results, is_empty, sense, objconstant, set_constant, constr_constant, rowranges, idxmap)
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
    # @show(user_settings)
end

hasresults(optimizer::Optimizer) = optimizer.hasresults

# MG: function to reset otimizer
function MOI.empty!(optimizer::Optimizer)
    COSMO.empty!(optimizer.inner)
    optimizer.hasresults = false
    optimizer.results = COSMO.Result{Float64}()
    optimizer.is_empty = true
    optimizer.sense = MOI.MinSense # model parameter, so needs to be reset
    optimizer.objconstant = 0.
    optimizer.set_constant = Float64[]
    optimizer.constr_constant = Float64[]   # the 5 in (3 * x1 + 2 * x2 + 5 ≤ 10 )
    optimizer.idxmap = MOIU.IndexMap()
    optimizer.rowranges = Dict{Int, UnitRange{Int}}()
    optimizer
end


# MG: check if optimizer object is empty
MOI.is_empty(optimizer::Optimizer) = optimizer.is_empty


struct UnsupportedObjectiveError <: Exception end


##############################
# MODEL --> SOLVER FORMAT FUNCTIONS
##############################

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; copy_names = false)
    copy_names && error("Copying names is not supported.")
    MOI.empty!(dest)
    idxmap = MOIU.IndexMap(dest, src)
    assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
    A,b, constr_constant, convexSets = processconstraints(dest, src, idxmap, dest.rowranges)
    COSMO.set!(dest.inner, P, q, A, b, convexSets, dest.inner.settings)
    dest.is_empty = false
    dest.idxmap = idxmap
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

# taken from https://github.com/JuliaOpt/SCS.jl

# Vectorized length for matrix dimension n
sympackedlen(n) = div(n*(n+1), 2)
# Matrix dimension for vectorized length n
sympackeddim(n) = div(isqrt(1+8n) - 1, 2)
trimap(i::Integer, j::Integer) = i < j ? trimap(j, i) : div((i-1)*i, 2) + j
trimapL(i::Integer, j::Integer, n::Integer) = i < j ? trimapL(j, i, n) : i + div((2n-j) * (j-1), 2)
function _sympackedto(x, n, mapfrom, mapto)
    @assert length(x) == sympackedlen(n)
    y = similar(x)
    for i in 1:n, j in 1:i
        y[mapto(i, j)] = x[mapfrom(i, j)]
    end
    y
end
sympackedLtoU(x, n=sympackeddim(length(x))) = _sympackedto(x, n, (i, j) -> trimapL(i, j, n), trimap)
sympackedUtoL(x, n=sympackeddim(length(x))) = _sympackedto(x, n, trimap, (i, j) -> trimapL(i, j, n))

function sympackedUtoLidx(x::AbstractVector{<:Integer}, n)
    y = similar(x)
    map = sympackedLtoU(1:sympackedlen(n), n)
    for i in eachindex(y)
        y[i] = map[x[i]]
    end
    y
end

orderval(val, s) = val
function orderval(val, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtoL(val, s.side_dimension)
end
orderidx(idx, s) = idx
function orderidx(idx, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtoLidx(idx, s.side_dimension)
end

# Scale coefficients depending on rows index
# rows: List of row indices
# coef: List of corresponding coefficients
# minus: if true, multiply the result by -1
# d: dimension of set
# rev: if true, we unscale instead (e.g. divide by √2 instead of multiply for PSD cone)
_scalecoef(rows, coef, minus, ::Type{<:MOI.AbstractSet}, d, rev) = minus ? -coef : coef
_scalecoef(rows, coef, minus, ::Union{Type{<:MOI.LessThan}, Type{<:MOI.Nonpositives}, Type{<:MOI.EqualTo}}, d, rev) = minus ? coef : -coef
function _scalecoef(rows, coef, minus, ::Type{MOI.PositiveSemidefiniteConeTriangle}, d, rev)
    scaling = minus ? -1 : 1
    scaling2 = rev ? scaling / √2 : scaling * √2
    output = copy(coef)
    diagidx = BitSet()
    for i in 1:d
        # vector index of diagonal entries of original matrix
        push!(diagidx, trimap(i, i))
    end
    # normalize rows
    rows = rows .- first(rows) .+ 1
    for i in 1:length(output)
        if rows[i] in diagidx
            output[i] *= scaling
        else
            output[i] *= scaling2
        end
    end
    output
end
# Unscale the coefficients in `coef` with respective rows in `rows` for a set `s` and multiply by `-1` if `minus` is `true`.
scalecoef(rows, coef, minus, s) = _scalecoef(rows, coef, minus, typeof(s), MOI.dimension(s), false)
# Unscale the coefficients in `coef` with respective rows in `rows` for a set of type `S` with dimension `d`
unscalecoef(rows, coef, S::Type{<:MOI.AbstractSet}, d) = _scalecoef(rows, coef, false, S, d, true)
#unscalecoef(rows, coef, S::Type{MOI.PositiveSemidefiniteConeTriangle}, d) = _scalecoef(rows, coef, false, S, sympackeddim(d), true)

function processconstraints(optimizer, src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}})

    m = mapreduce(length, +, values(rowranges), init=0)

    b = zeros(Float64, m)
    constant = zeros(Float64, m)
    I = Int[]
    J = Int[]
    V = Float64[]
    COSMOconvexSets = Array{COSMO.AbstractConvexSet{Float64}}(undef, 0)
    set_constant = zeros(Float64, m)
    # loop over constraints and modify A, l, u and constants
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        processconstraints!((I, J, V), b, COSMOconvexSets, constant, set_constant, src, idxmap, rowranges, F, S)
    end
    optimizer.set_constant = set_constant
    # subtract constant from right hand side
    n = MOI.get(src, MOI.NumberOfVariables())
    A = sparse(I, J, V, m, n)
    return A, b, constant, COSMOconvexSets
end

function processconstraints!(triplets::SparseTriplets, b::Vector, COSMOconvexSets::Array{COSMO.AbstractConvexSet{Float64}}, constant::Vector{Float64}, set_constant::Vector{Float64},
        src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}},
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    # loop over all constraints of same (F,S)
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        # handle cases where the sign of rows of A an b has to be flipped
        # isa(s, Union{MOI.GreaterThan, MOI.Nonnegatives, SOC, PSD}) ? FLIP_SIGN = true : FLIP_SIGN = false

        rows = constraint_rows(rowranges, idxmap[ci])
        processConstant!(b, rows, f, s)
        constant[rows] = b[rows]
        if typeof(s) <: MOI.AbstractScalarSet
            set_constant[rows] =  MOIU.getconstant(s)
        end

        processConstraint!(triplets, f, rows, idxmap, s)
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
function processConstant!(b, row::Int, f::AffineConvertible, s)
    b[row] = scalecoef(row, constant(f), false, s)
    nothing
end

function processConstant!(b, rows::UnitRange{Int}, f::VectorOfVariables, s)
    b[rows] .= 0
    nothing
end

# process constant for functions VectorAffineFunction{Float64}
function processConstant!(b, rows::UnitRange{Int}, f::VectorAffine, s)
    for (i, row) in enumerate(rows)
        b[row] = scalecoef(row, f.constants[i], false, s)
    end
    nothing
end

# process function like f(x)= x1
function processConstraint!(triplets::SparseTriplets, f::MOI.SingleVariable, row::Int, idxmap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    col = idxmap[f.variable].value
    push!(I, row)
    push!(J, col)
    push!(V, scalecoef(row, 1, true, s))
    nothing
end

# process function like f(x) = coeff*x_1
function processConstraint!(triplets::SparseTriplets, f::MOI.ScalarAffineFunction, row::Int, idxmap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    for term in f.terms
        var = term.variable_index
        coeff = term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, scalecoef(row, coeff, true, s))
    end
end

function processConstraint!(triplets::SparseTriplets, f::MOI.VectorOfVariables, rows::UnitRange{Int}, idxmap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    cols = zeros(length(rows))
    for (i, var) in enumerate(f.variables)
        cols[i] = idxmap[var].value
    end
    row_start = first(rows) - 1
    # @show(row_start, orderidx(collect(1:length(rows)), s), cols, scalecoef(rows, orderval(ones(length(rows)), s), true, s))
    append!(I, rows)#row_start .+ orderidx(collect(1:length(rows)), s))
    append!(J, cols)
    append!(V, scalecoef(rows, orderval(ones(length(rows)), s), true, s))
    nothing
end

# FIXME: This won't work with the loop because of the way scalecoef is written
# process function like f(x) = a'*x where x[1:n]
function processConstraint!(triplets::SparseTriplets, f::MOI.VectorAffineFunction, rows::UnitRange{Int}, idxmap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    for term in f.terms
        row = rows[term.output_index]
        var = term.scalar_term.variable_index
        coeff = term.scalar_term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, scalecoef(rows, coeff, true, s))
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

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.Nonpositives)
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
# GET / SET : OPTIMIZER ATTRIBUTES
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
    if optimizer.sense == MOI.MaxSense
        return -rawobj - optimizer.objconstant
    else
        return rawobj + optimizer.objconstant
    end
end

# Get solve time
MOI.get(optimizer::Optimizer, a::MOI.SolveTime) = optimizer.results.times.solverTime

MOI.get(optimizer::Optimizer, a::MOI.SolverName) = "COSMO"

# Get Termination Status
function MOI.get(optimizer::Optimizer, a::MOI.TerminationStatus)
    # Note that the :Dual_infeasible and :Primal_infeasible are mapped to MOI.Success
    # because COSMO can return a proof of infeasibility. For the same reason,
    # :Primal_infeasible_inaccurate is mapped to MOI.AlmostSuccess
    status = optimizer.results.status
    if status == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif status == :Dual_infeasible
        MOI.Success
    elseif status == :Primal_infeasible
        MOI.Success
    elseif status == :Max_iter_reached
        MOI.IterationLimit
    elseif status == :Solved
        MOI.Success
    end
end

# Get Primal Status
function MOI.get(optimizer::Optimizer, a::MOI.PrimalStatus)
    status = optimizer.results.status
    if status == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif status == :Primal_infeasible
        MOI.InfeasiblePoint
    elseif status == :Solved
        MOI.FeasiblePoint
    elseif status == :Dual_infeasible
        MOI.InfeasibilityCertificate
    else
        MOI.UnknownResultStatus
    end
end

# Get Dual Status
function MOI.get(optimizer::Optimizer, a::MOI.DualStatus)
    status = optimizer.results.status
    if status == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif status == :Dual_infeasible
        MOI.InfeasiblePoint
    elseif status == :Primal_infeasible
        MOI.InfeasibilityCertificate
    elseif status == :Solved
        MOI.FeasiblePoint
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate, :Non_convex (TODO: good idea? use COSMO.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end

# function constraint_primal(optimizer::Optimizer, rows, S::Type{<:MOI.AbstractSet})
#      @show(optimizer.inner.p.b, optimizer.constr_constant, optimizer.results.s )
#      coeff = optimizer.inner.p.b[rows] + optimizer.constr_constant[rows] - optimizer.results.s[rows]
#        return unscalecoef(rows, reorderval(optimizer.sol.slack[offset .+ rows], S), S, length(rows)), S)

# end

_unshift(optimizer::Optimizer, offset, value, s) = value
_unshift(optimizer::Optimizer, offset, value, s::Type{<:MOI.AbstractScalarSet}) = value + optimizer.set_constant[offset]

function MOI.get(optimizer::Optimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    # offset = constroffset(optimizer, ci)
    # rows = constrrows(optimizer, ci)
    # _unshift(optimizer, offset, unscalecoef(rows, reorderval(optimizer.sol.slack[offset .+ rows], S), S, length(rows)), S)
    rows = constraint_rows(optimizer.rowranges, optimizer.idxmap[ci])
    c_primal = unscalecoef(rows, optimizer.results.s[rows], S, length(rows))
    # (typeof(c_primal) <: AbstractArray && length(c_primal) == 1) && (c_primal = first(c_primal))
    return _unshift(optimizer, rows, c_primal, S)
end

function MOI.get(optimizer::Optimizer, ::MOI.ConstraintDual, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    rows = constraint_rows(optimizer.rowranges, optimizer.idxmap[ci])
    return unscalecoef(rows, optimizer.results.y[rows], S, length(rows))
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
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Union{VectorOfVariables, VectorAffine}}, ::Type{<:SupportedVectorSets}) = true
# MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Union{VectorOfVariables,VectorAffine}}, ::Type{<:PSD}) = true



