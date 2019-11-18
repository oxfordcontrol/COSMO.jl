# This interface was partially derived from the OSQP MathOptInterface written by T. Koolen and B. Stellato
# https://github.com/oxfordcontrol/OSQP.jl/blob/master/src/MathOptInterfaceOSQP.jl
# certain utility function are taken from the SCS  MathOptInterface:
# https://github.com/JuliaOpt/SCS.jl
using MathOptInterface

const MOI = MathOptInterface
const MOIU = MOI.Utilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const SparseTriplets = Tuple{Vector{Int}, Vector{Int}, Vector{<:Any}}

const SingleVariable = MOI.SingleVariable
const Affine = MOI.ScalarAffineFunction{Float64}
const Quadratic = MOI.ScalarQuadraticFunction{Float64}
const VectorOfVariables = MOI.VectorOfVariables
const VectorAffine = MOI.VectorAffineFunction{Float64}

const Interval = MOI.Interval{Float64}
const LessThan = MOI.LessThan{Float64}
const GreaterThan = MOI.GreaterThan{Float64}
const EqualTo = MOI.EqualTo{Float64}
const IntervalConvertible = Union{LessThan, GreaterThan, EqualTo, Interval}


const Zeros = MOI.Zeros
#Nonnegatives = MOI.Nonnegatives
const Nonpositives = MOI.Nonpositives
const SOC = MOI.SecondOrderCone
const PSDS = MOI.PositiveSemidefiniteConeSquare
const PSDT = MOI.PositiveSemidefiniteConeTriangle
const PSD = Union{MOI.PositiveSemidefiniteConeSquare,MOI.PositiveSemidefiniteConeTriangle}
const SupportedVectorSets = Union{Zeros, MOI.Nonnegatives, Nonpositives, SOC, PSDS, PSDT, MOI.ExponentialCone, MOI.DualExponentialCone, MOI.PowerCone, MOI.DualPowerCone}

#export sortSets, assign_constraint_row_ranges!, processconstraints, constraint_rows, processobjective, processlinearterms!, symmetrize!, processconstraints!, constant, processconstant!, processlinearpart!, processconstraintset!
export Optimizer

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
        set_MOI_default_settings!(inner.settings)
        length(user_settings) > 0 && set_user_settings!(inner, user_settings)
        hasresults = false
        results = COSMO.Result{Float64}()
        is_empty = true
        sense = MOI.MIN_SENSE
        objconstant = 0.
        set_constant = Float64[]
        constr_constant = Float64[]
        rowranges = Dict{Int, UnitRange{Int}}()
        idxmap = MOIU.IndexMap()
        new(inner, hasresults, results, is_empty, sense, objconstant, set_constant, constr_constant, rowranges, idxmap)
    end
end

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
function MOI.empty!(optimizer::Optimizer)
    COSMO.empty_model!(optimizer.inner)
    optimizer.hasresults = false
    optimizer.results = COSMO.Result{Float64}()
    optimizer.is_empty = true
    optimizer.sense = MOI.MIN_SENSE # model parameter, so needs to be reset
    optimizer.objconstant = 0.
    optimizer.set_constant = Float64[]
    optimizer.constr_constant = Float64[]   # the 5 in (3 * x1 + 2 * x2 + 5 ≤ 10 )
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
    assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    processconstraints!(dest, src, idxmap, dest.rowranges)
    pre_allocate_variables!(dest.inner)

    dest.is_empty = false
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
    # LOCs = MOI.get(src, MOI.ListOfConstraints())
    # sort!(LOCs, by=x-> sort_sets(x[2]))
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
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
sort_sets(s::Union{Type{<: MOI.EqualTo}, Type{<: MOI.Zeros}}) = 1
sort_sets(s::Union{Type{ <: LessThan}, Type{ <: GreaterThan}, Type{ <: MOI.Nonnegatives}, Type{ <: MOI.Nonpositives}}) = 2
sort_sets(s::Type{MOI.Interval}) = 3
sort_sets(s::Type{<: MOI.AbstractSet}) = 4

# this has to be in the right order
# returns a Dictionary that maps a constraint index to a range
function assign_constraint_row_ranges!(rowranges::Dict{Int, UnitRange{Int}}, idxmap::MOIU.IndexMap, src::MOI.ModelLike)
    startrow = 1
    LOCs = MOI.get(src, MOI.ListOfConstraints())
    sort!(LOCs, by=x-> sort_sets(x[2]))
    for (F, S) in LOCs
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

# this is handled explicitly in `copy_to` before the objective function is processed
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveSense, sense) = nothing

function allocate_cost_variables!(optimizer::Optimizer)
    m, n = size(optimizer.inner.p.A)
    optimizer.inner.p.P = spzeros(n, n)
    optimizer.inner.p.q = zeros(n)
    return nothing
end

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{<:Union{MOI.SingleVariable, Affine, Quadratic}}) = true
function MOIU.load(dest::Optimizer, ::MOI.ObjectiveFunction, f::MOI.SingleVariable)
    idxmap = dest.idxmap
    n = MOI.get(dest, MOI.NumberOfVariables())
    dest.inner.p.q = zeros(n)
    dest.inner.p.q[idxmap[f.variable].value] = 1
    dest.inner.p.P = spzeros(n, n)
    dest.objconstant = 0.
    dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
    return nothing
end

function MOIU.load(dest::Optimizer, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction{Float64})
    idxmap = dest.idxmap
    n = MOI.get(dest, MOI.NumberOfVariables())
    dest.inner.p.P = spzeros(n, n)

    dest.inner.p.q = zeros(n)
    processlinearterms!(dest.inner.p.q, f.terms)# idxmap)

    dest.objconstant = f.constant
    dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
    return nothing
end

function MOIU.load(dest::Optimizer, ::MOI.ObjectiveFunction, f::MOI.ScalarQuadraticFunction{Float64})
    idxmap = dest.idxmap
    n = MOI.get(dest, MOI.NumberOfVariables())

    I = [Int(idxmap[term.variable_index_1].value) for term in f.quadratic_terms]
    J = [Int(idxmap[term.variable_index_2].value) for term in f.quadratic_terms]
    V = [term.coefficient for term in f.quadratic_terms]
    symmetrize!(I, J, V)
    dest.inner.p.P = sparse(I, J, V, n, n)

    dest.inner.p.q = zeros(n)
    processlinearterms!(dest.inner.p.q, f.affine_terms)#, idxmap)
    dest.objconstant = f.constant

    dest.objconstant = apply_sense!(dest.sense, dest.inner.p.P, dest.inner.p.q, dest.objconstant)
    return nothing
end

function MOIU.load(dest::Optimizer, ::MOI.ObjectiveFunction, f)
    throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction{typeof(f)}()))
end

function apply_sense!(sense::MOI.OptimizationSense, P::AbstractMatrix, q::AbstractVector, c::Float64)
    if sense == MOI.MAX_SENSE
        LinearAlgebra.rmul!(P, -1)
        LinearAlgebra.rmul!(q, -1)
        c = -c
    end
    return c
end

function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm}, idxmapfun::Function = identity)
    for term in terms
        var = term.variable_index
        coeff = term.coefficient
        q[var.value] += coeff
    end
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

# The following utility functions are from https://github.com/JuliaOpt/SCS.jl

# Scale coefficients depending on rows index
# rows: List of row indices
# coef: List of corresponding coefficients
# minus: if true, multiply the result by -1
# d: dimension of set
# rev: if true, we unscale instead (e.g. divide by √2 instead of multiply for PSD cone)
_scalecoef(rows, coef, minus, ::Type{<:MOI.AbstractSet}, rev) = minus ? -coef : coef
_scalecoef(rows, coef, minus, ::Union{Type{<:MOI.LessThan}, Type{<:MOI.Nonpositives}, Type{<:MOI.EqualTo}}, rev) = minus ? coef : -coef
function _scalecoef(rows, coef, minus, ::Type{MOI.PositiveSemidefiniteConeTriangle}, rev)
    scaling = minus ? -1 : 1
    scaling2 = rev ? scaling / √2 : scaling * √2
    output = copy(coef)
    idx = 0
    for i in 1:length(output)
        # See https://en.wikipedia.org/wiki/Triangular_number#Triangular_roots_and_tests_for_triangular_numbers
        is_diagonal_index = isinteger(sqrt(8*rows[i] + 1))
        if is_diagonal_index
            output[i] *= scaling
        else
            output[i] *= scaling2
        end
    end
    output
end
# Unscale the coefficients in `coef` with respective rows in `rows` for a set `s` and multiply by `-1` if `minus` is `true`.
scalecoef(rows, coef, minus, s) = _scalecoef(rows, coef, minus, typeof(s), false)
# Unscale the coefficients in `coef` with respective rows in `rows` for a set of type `S`
unscalecoef(rows, coef, S::Type{<:MOI.AbstractSet}) = _scalecoef(rows, coef, false, S, true)

# This helper function is to provide scale- and unscalecoef with the nomimal rows in the case of at PSDTriangle constraint
# this is because for _scalecoef the actual row in A doesn't matter, what matters is the row in the upper triangle matrix
nom_rows(rows, s::Type{MOI.PositiveSemidefiniteConeTriangle}) = 1:length(rows)
nom_rows(rows, s::Type{<:MOI.AbstractSet}) = rows

output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)

function processconstraints!(optimizer::Optimizer, src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}})

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

    # in a second step: merge some of the sets
    convex_sets = COSMO.merge_sets(COSMOconvexSets)

    # store constraint matrices in inner model
    model = optimizer.inner
    model.p.A = A
    model.p.b = b
    model.p.C = CompositeConvexSet(deepcopy(convex_sets))
    model.p.model_size = [size(A, 1); size(A,2 )]
    optimizer.constr_constant = constant
    return nothing
end

function processconstraints!(triplets::SparseTriplets, b::Vector, COSMOconvexSets::Array{COSMO.AbstractConvexSet{Float64}}, constant::Vector{Float64}, set_constant::Vector{Float64},
        src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}},
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    # loop over all constraints of same (F,S)
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        rows = constraint_rows(rowranges, idxmap[ci])
        processConstant!(b, rows, f, s)
        constant[rows] = b[rows]
        if typeof(s) <: MOI.AbstractScalarSet && !(typeof(s) <: MOI.Interval)
            set_constant[rows] =  MOI.constant(s)
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
function processConstant!(b, row::Int, f::Affine, s)
    b[row] = scalecoef(row, constant(f), false, s)
    nothing
end

function processConstant!(b, rows::UnitRange{Int}, f::VectorOfVariables, s)
    b[rows] .= 0
    nothing
end

# process constant for functions VectorAffineFunction{Float64}
function processConstant!(b, rows::UnitRange{Int}, f::VectorAffine, s)
        b[rows] = scalecoef(nom_rows(rows, typeof(s)), f.constants, false, s)
    nothing
end


# process function like f(x) = coeff*x_1
function processConstraint!(triplets::SparseTriplets, f::MOI.ScalarAffineFunction, row::Int, idxmap::MOIU.IndexMap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    for term in f.terms
        push!(I, row)
        push!(J, idxmap[term.variable_index].value)
        push!(V, scalecoef(row, term.coefficient, true, s))
    end
end

function processConstraint!(triplets::SparseTriplets, f::MOI.VectorOfVariables, rows::UnitRange{Int}, idxmap::MOIU.IndexMap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    cols = zeros(length(rows))
    for (i, var) in enumerate(f.variables)
        cols[i] = idxmap[var].value
    end
    append!(I, rows)
    append!(J, cols)
    append!(V, scalecoef(nom_rows(rows, typeof(s)), ones(length(rows)), true, s))
    nothing
end

function processConstraint!(triplets::SparseTriplets, f::MOI.VectorAffineFunction, rows::UnitRange{Int}, idxmap::MOIU.IndexMap, s::MOI.AbstractSet)
    (I, J, V) = triplets
    A = sparse(output_index.(f.terms), variable_index_value.(f.terms), coefficient.(f.terms))
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    A_I, A_J, A_V = findnz(A)

    offset = first(rows) - 1

    append!(J, A_J)
    append!(V, scalecoef(A_I, A_V, true, s))
    append!(I, A_I .+ offset)
end


##############################
# PROCESS SETS
##############################

# process the following sets Union{LessThan, GreaterThan, EqualTo}
function processSet!(b::Vector, row::Int, cs, s::LessThan)
    b[row] += s.upper
    push!(cs, COSMO.Nonnegatives{Float64}(1))
    nothing
end
function processSet!(b::Vector, row::Int, cs, s::GreaterThan)
    b[row] -= s.lower
    push!(cs, COSMO.Nonnegatives{Float64}(1))
    nothing
end
function processSet!(b::Vector, row::Int, cs, s::EqualTo)
    b[row] += s.value
    push!(cs, COSMO.ZeroSet{Float64}(1))
    nothing
end

function processSet!(b::Vector, row::Int, cs, s::MOI.Interval)
    push!(cs, COSMO.Box{Float64}([s.lower], [s.upper]))
    nothing
end

# process the following sets Zeros
function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::Zeros)
    push!(cs, COSMO.ZeroSet{Float64}(length(rows)))
    nothing
end

# process the following sets Union{Zeros, Nonnegatives, Nonpositives}
function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.Nonnegatives)
    push!(cs, COSMO.Nonnegatives{Float64}(length(rows)))
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.Nonpositives)
    push!(cs, COSMO.Nonnegatives{Float64}(length(rows)))
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::SOC)
    push!(cs, COSMO.SecondOrderCone{Float64}(length(rows)))
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.PositiveSemidefiniteConeSquare)
    push!(cs, COSMO.PsdCone{Float64}(length(rows)))
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.PositiveSemidefiniteConeTriangle)
    push!(cs, COSMO.PsdConeTriangle{Float64}(length(rows)))
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.ExponentialCone)
    push!(cs, COSMO.ExponentialCone{Float64}())
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.DualExponentialCone)
    push!(cs, COSMO.DualExponentialCone{Float64}())
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.PowerCone)
    push!(cs, COSMO.PowerCone{Float64}(s.exponent))
    nothing
end

function processSet!(b::Vector, rows::UnitRange{Int}, cs, s::MOI.DualPowerCone)
    push!(cs, COSMO.DualPowerCone{Float64}(s.exponent))
    nothing
end


# to reduce function calls combine, individual ZeroSets, Nonnegatives and Box constraints
# This assumes that C is ordered by set type
function merge_sets(C::Array{COSMO.AbstractConvexSet{Float64}})

    z_ind = findall(set -> typeof(set) <: COSMO.ZeroSet, C)
    nn_ind = findall(set -> typeof(set) <: COSMO.Nonnegatives, C)
    box_ind = findall(set -> typeof(set) <: COSMO.Box, C)
    other_ind = findall( x-> !in(x, union(z_ind, nn_ind, box_ind)), collect(1:1:length(C)))

    # if none of these sets are present, do nothing
    length(other_ind) == length(C) && return C

    length(z_ind) > 0 ? (num_z = 1) : (num_z = 0)
    length(nn_ind) > 0 ? (num_nn = 1) : (num_nn = 0)
    length(box_ind) > 0 ? (num_box = 1) : (num_box = 0)


    num_merged_sets = length(other_ind) + num_z + num_nn + num_box
    merged_sets = Array{COSMO.AbstractConvexSet{Float64}}(undef, num_merged_sets)
    set_ind = 1
    length(z_ind) > 0 && (set_ind = COSMO.merge_set!(merged_sets, z_ind, C, set_ind, COSMO.ZeroSet{Float64}))
    length(nn_ind) > 0 && (set_ind = COSMO.merge_set!(merged_sets, nn_ind, C, set_ind, COSMO.Nonnegatives{Float64}))
    length(box_ind) > 0 && (set_ind = COSMO.merge_box!(merged_sets, box_ind, C, set_ind))

    for other in other_ind
        merged_sets[set_ind] = C[other]
        set_ind += 1
    end
    return merged_sets
end


function merge_set!(merged_sets::Array{COSMO.AbstractConvexSet{Float64}, 1}, ind::Array{Int64, 1}, C::Array{<: COSMO.AbstractConvexSet, 1}, set_ind::Int64, set_type::DataType)
        if length(ind) > 1
            combined_dim = sum(x -> x.dim, C[ind])
        else
            combined_dim = C[ind[1]].dim
        end
        merged_sets[set_ind] = set_type(combined_dim)
        return set_ind + 1
end

function merge_box!(merged_sets::Array{COSMO.AbstractConvexSet{Float64}, 1}, ind::Array{Int64, 1}, C::Array{<: COSMO.AbstractConvexSet{Float64}, 1}, set_ind::Int64)
        if length(ind) > 1
            combined_dim = sum(x -> x.dim, C[ind])
            l = zeros(Float64, combined_dim)
            u = zeros(Float64, combined_dim)
            row = 1
            for box in C[ind]
                l[row:row + box.dim - 1] = box.l
                u[row:row + box.dim - 1] = box.u
                row += box.dim
            end
            merged_sets[set_ind] = COSMO.Box(l, u)
        else
            merged_sets[set_ind] = C[ind[1]]
        end

        return set_ind + 1
end



function pass_attributes!(dest::Optimizer, src::MOI.ModelLike, idxmap::MOIU.IndexMap, pass_attr::Function=MOI.set)
    copy_names = false

    # Copy model attributes, e.g. ObjectiveFunction
    attrs = MOI.get(src, MOI.ListOfModelAttributesSet())
    MOIU._pass_attributes(dest, src, copy_names, idxmap, attrs, tuple(), tuple(), tuple(), MOIU.load)

    # Copy variable attributes, e.g. VariablePrimalStart
    var_attr = MOI.get(src, MOI.ListOfVariableAttributesSet())
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    vis_dest = map(vi -> idxmap[vi], vis_src)
    MOIU._pass_attributes(dest, src, copy_names, idxmap, var_attr, (VI,), (vis_src,), (vis_dest,), MOIU.load)

    # Copy constraint attributes, e.g. ConstraintPrimalStart
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        attrs = MOI.get(src, MOI.ListOfConstraintAttributesSet{F, S}())
        cis_dest = map(ci -> idxmap[ci], cis_src)
        MOIU._pass_attributes(dest, src, copy_names, idxmap, attrs, (CI{F, S},), (cis_src,), (cis_dest,), MOIU.load)
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

MOI.get(optimizer::Optimizer, a::MOI.SolveTime) = optimizer.results.times.solver_time
MOI.get(optimizer::Optimizer, a::MOI.SolverName) = "COSMO"

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
    if !(1 <= a.N <= MOI.get(optimizer, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

    status = optimizer.results.status
    if status == :Primal_infeasible
        return MOI.INFEASIBLE_POINT
    elseif status == :Solved
        return MOI.FEASIBLE_POINT
    elseif status == :Dual_infeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        return MOI.NO_SOLUTION
    end
end

# Get Dual Status
function MOI.get(optimizer::Optimizer, a::MOI.DualStatus)
    if !(1 <= a.N <= MOI.get(optimizer, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

    status = optimizer.results.status
    if status == :Dual_infeasible
        return MOI.INFEASIBLE_POINT
    elseif status == :Primal_infeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif status == :Solved
        return MOI.FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    return String(optimizer.results.status)
end




_unshift(optimizer::Optimizer, offset, value, s) = value
_unshift(optimizer::Optimizer, offset, value, s::Type{<:MOI.AbstractScalarSet}) = value + optimizer.set_constant[offset]
_shift(optimizer::Optimizer, offset, value, s) = value
_shift(optimizer::Optimizer, offset, value, s::Type{<:MOI.AbstractScalarSet}) = value - optimizer.set_constant[offset]

function MOI.get(optimizer::Optimizer, a::MOI.ConstraintPrimal, ci_src::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, a)
    # offset = constroffset(optimizer, ci)
    # rows = constrrows(optimizer, ci)
    rows = COSMO.constraint_rows(optimizer.rowranges, ci_src)
    c_primal = unscalecoef(nom_rows(rows, S), optimizer.results.s[rows], S)
    # (typeof(c_primal) <: AbstractArray && length(c_primal) == 1) && (c_primal = first(c_primal))
    return _unshift(optimizer, rows, c_primal, S)
end

function MOI.get(optimizer::Optimizer, a::MOI.ConstraintDual, ci_src::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, a)
    rows = constraint_rows(optimizer.rowranges, ci_src)
    return unscalecoef(nom_rows(rows, S), optimizer.results.y[rows], S)
end

## Variable attributes:
function MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::VI)
    MOI.check_result_index_bounds(optimizer, a)
    return optimizer.results.x[vi.value]
end

## Warm starting:
MOI.supports(::Optimizer, a::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}) = true
function MOIU.load(optimizer::Optimizer, a::MOI.VariablePrimalStart, vi::VI, value::Real)
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    COSMO.warm_start_primal!(optimizer.inner, value, vi.value)
end
MOIU.load(optimizer::Optimizer, a::MOI.VariablePrimalStart, vi::VI, value::Nothing) = nothing

MOI.supports(::Optimizer, a::MOI.ConstraintPrimalStart, ::Type{<:MOI.ConstraintIndex}) = true
function MOIU.load(optimizer::Optimizer, a::MOI.ConstraintPrimalStart, ci::CI{<:MOI.AbstractFunction, S}, value) where S <: MOI.AbstractSet
    (value == nothing || isa(value, Array{Nothing, 1})) && return nothing
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    rows = constraint_rows(optimizer.rowranges, ci)
    # this undoes the scaling and shifting that was used in get(MOI.ConstraintPrimal)
    # Off-diagonal entries of slack variable of a PSDTriangle constraint has to be scaled by sqrt(2)
    value = _shift(optimizer, rows, value, S)
    COSMO.warm_start_slack!(optimizer.inner, _scalecoef(nom_rows(rows, S), value, false, S, false), rows)
end

MOI.supports(::Optimizer, a::MOI.ConstraintDualStart, ::Type{<:MOI.ConstraintIndex}) = true
function MOIU.load(optimizer::Optimizer, a::MOI.ConstraintDualStart, ci::CI{<:MOI.AbstractFunction, S}, value) where S <: MOI.AbstractSet
    (value == nothing || isa(value, Array{Nothing, 1})) && return nothing
    MOI.is_empty(optimizer) && throw(MOI.CannotSetAttribute(a))
    rows = constraint_rows(optimizer.rowranges, ci)
    # this undoes the scaling that was used in get(MOI.ConstraintDual)
    # Off-diagonal entries of dual variable of a PSDTriangle constraint has to be scaled by sfqrt(2)
    COSMO.warm_start_dual!(optimizer.inner, _scalecoef(nom_rows(rows, S), value, false, S, false), rows)
    nothing
end


 MOI.supports(::Optimizer, ::MOI.RawParameter) = true
function MOI.set(optimizer::Optimizer, p::MOI.RawParameter, value)
    setfield!(optimizer.inner.settings, Symbol(p.name), value)
end

function MOI.get(optimizer::Optimizer, p::MOI.RawParameter)
    if in(Symbol(p.name), fieldnames(typeof(optimizer.inner.settings)))
        return getfield(optimizer.inner.settings, Symbol(p.name))
    end
    error("RawParameter with name $(p.name) is not set.")
end


# MOI.Silent == COSMO.Settings.verbose = false
MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.inner.settings.verbose = !value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = !optimizer.inner.settings.verbose

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(optimizer::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(optimizer, MOI.RawParameter("time_limit"), Float64(value))
end

function MOI.set(optimizer::Optimizer, attr::MOI.TimeLimitSec, ::Nothing)
    MOI.set(optimizer, MOI.RawParameter("time_limit"), 0.0)
end
function MOI.get(optimizer::Optimizer, ::MOI.TimeLimitSec)
    return MOI.get(optimizer, MOI.RawParameter("time_limit"))
end

##############################
# SUPPORT OBJECTIVE AND SET FUNCTIONS
##############################


## Supported Constraints
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Affine}, ::Type{<:IntervalConvertible}) = true
MOI.supports_constraint(optimizer::Optimizer, ::Type{<:Union{VectorOfVariables, VectorAffine}}, ::Type{<:SupportedVectorSets}) = true
