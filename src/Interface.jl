# User facing functions / structs:

# ------------------------------
# Structs: (defined in Types.jl)
# ------------------------------
# 1. QOCS.Model
# 2. QOCS.Settings
# 3. QOCS.Results
# 4. QOCS.ConvexSet (Abstract)
# 4a. Subtypes of ConvexSet: ZeroCone, NonNegativeOrthant, Interval (Box), SecondOrderCone, PositiveSemidefiniteCone <: QOCS.ConvexSet, User defined ones
# Each Subtype of ConvexSet needs an associated projection function project!(set::Type{<:ConvexSet}). We'll provide the standard ones, and a isInDualCone->Bool, isInRecessionCone->Bool

# ------------------------------
# Functions: (defined in Interface.jl)
# ------------------------------

# 1. assemble!() --> To copy problem data and warm start data into the model
# Idea: two assemble! functions:
# 1a. Low level: Pass the problem data exactly in the current problem format(A,q,A,b,K) with K sedumi style
# 1b. High Level: Pass P,q directly, pass constraints as an array of ConvexSet Subtypes

# 2. optimize!(model::QOCS.Model,settings::QOCS.settings) --> Calls the solver based on the model data and returns a QOCS.Results object
# 3. reset!() --> Empty functions to reset each struct QOCS.Model, QOCS.Settings, QOCS.Results objects

# ------------------------------
# Open questions:
# ------------------------------
# 1. How best handle warm starting?
# 2. How best handle update of problem data without re-assambling the whole model?
# 3. Ideally work directly on the QOCS.Model data without copying anything (that means that matrixes will be changed through scaling) --> Problem?

# function assemble!(model::QOCS.Model,P,q,A,b,K)
#   model.P = P
#   model.q = q
#   model.A = A
#   model.b = b
#   model.K = K
#   model.m, model.n = size(A)
# end

function assemble!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::AbstractVector{<:Real},constraints::Array{QOCS.Constraint},x0::Union{Vector{Float64}, Nothing} = nothing, y0::Union{Vector{Float64}, Nothing} = nothing)
  # convert inputs
  P[:,:] = convert(SparseMatrixCSC{Float64,Int64},P)
  q[:] = convert(Vector{Float64},q)

  # model.Flags.INFEASIBILITY_CHECKS = checkConstraintFunctions(constraints)

  mergeConstraints!(constraints)
  model.P = P
  model.q = q
  n = length(q)
  m = sum(map(x->x.dim,map(x->x.convexSet,constraints)))

  model.m = m
  model.n = n
  model.A = spzeros(Float64,m,n)
  model.b = spzeros(Float64,m)

  model.x0 = zeros(Float64,n)
  model.y0 = zeros(Float64,m)

  # merge and sort the constraint sets
  sort!(constraints,by=x->sortSets(x.convexSet))
  rowNum = 1
  for con in constraints
    processConstraint!(model,rowNum,con.A,con.b,con.convexSet)
    rowNum += con.convexSet.dim
  end

  # save the convex sets inside the model
  model.convexSets = map(x->x.convexSet,constraints)

  # if user provided warm starting variables, update model
  warmStart!(model,x0 = x0,y0 = y0)

  nothing
end

function warmStart!(model::QOCS.Model; x0::Union{Vector{Float64}, Nothing} = nothing, y0::Union{Vector{Float64}, Nothing} = nothing)
    x0 isa Vector{Float64} && (model.x0 = x0)
    y0 isa Vector{Float64} && (model.y0 = y0)
end

function assemble!(model::QOCS.Model,P::Real,q::Real,constraints::Array{QOCS.Constraint})
  Pm = spzeros(1,1)
  qm = zeros(1)
  Pm[1,1] = convert(Float64,P)
  qm[1] = convert(Float64,q)
  assemble!(model,Pm,qm,constraints)
end

assemble!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::AbstractMatrix{<:Real},constraints::Array{QOCS.Constraint}) = assemble!(model,P,vec(q),constraints)
assemble!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::Real,constraints::Array{QOCS.Constraint}) = assemble!(model,P,[q],constraints)


# merge zeros sets and nonnegative sets
function mergeConstraints!(constraints::Array{QOCS.Constraint})
  # handle zeros sets
  ind = findall(set->typeof(set) == Zeros,map(x->x.convexSet,constraints))
  if length(ind) > 1
    M = mergeZeros(constraints[ind])
    deleteat!(constraints,ind)
    push!(constraints,M)
  end

  # handle nonnegative sets
  ind = findall(set->typeof(set) == Nonnegatives,map(x->x.convexSet,constraints))
  if length(ind) > 1
    M = mergeNonnegatives(constraints[ind])
    deleteat!(constraints,ind)
    push!(constraints,M)
  end
  nothing
end

function mergeZeros(constraints::Array{QOCS.Constraint})
  m = sum(x->x.dim,map(x->x.convexSet,constraints))
  n = size(constraints[1].A,2)
  A = spzeros(m,n)
  b = zeros(m)

  s = 1
  e = 0
  for cons in constraints
    e = s + cons.convexSet.dim -1
    A[s:e,:] = cons.A
    b[s:e,:] = cons.b
    s = e + 1
  end

  return M = QOCS.Constraint(A,b,Zeros())
end

function mergeNonnegatives(constraints::Array{QOCS.Constraint})
  m = sum(x->x.dim,map(x->x.convexSet,constraints))
  n = size(constraints[1].A,2)
  A = spzeros(m,n)
  b = zeros(m)

  s = 1
  e = 0
  for cons in constraints
    e = s + cons.convexSet.dim -1
    A[s:e,:] = cons.A
    b[s:e,:] = cons.b
    s = e + 1
  end

  return M = QOCS.Constraint(A,b,Nonnegatives())
end


function sortSets(C::AbstractConvexSet)
  C = typeof(C)
  (C <: Zeros) && return 1
  (C <: Nonnegatives) && return 2
  (C <: Box) && return 3
  (C <: SecondOrderCone) && return 4
  (C <: PositiveSemidefiniteCone) && return 5
  return 6
end

# transform A*x + b in {0}, to A*x + s == b, s in {0}
function processConstraint!(model::QOCS.Model,rowNum::Int64,A::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},b::AbstractVector{<:Real},C::AbstractConvexSet)
  s = rowNum
  e = rowNum + C.dim - 1
  model.A[s:e,:] = -A
  model.b[s:e,:] = b
  C.indices = s:e
end
