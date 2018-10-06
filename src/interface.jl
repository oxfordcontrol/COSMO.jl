"""
    assemble!(model,P,q,constraint(s),[x0,y0])

Assembles a `QOCS.Model` with a cost function defind by `P` and `q`, and a number of `constraints`.

The positive semidefinite matrix `P` and vector `q` are used to specify the cost function of the optimization problem:

```
min   1/2 x'Px + q'x
s.t.  Ax + b ∈ C
```
`constraints` is a `QOCS.Constraint` or an array of `QOCS.Constraint` objects that are used to describe the constraints on `x`.

---
The optinal arguments `x0` and `y0` can be used to provide the solver with warm starting values for the primal variable `x` and the dual variable `y`.
"""
function assemble!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::AbstractVector{<:Real},constraints::Union{QOCS.Constraint,Array{QOCS.Constraint,1}},x0::Union{Vector{Float64}, Nothing} = nothing, y0::Union{Vector{Float64}, Nothing} = nothing)
  # convert inputs
  P[:,:] = convert(SparseMatrixCSC{Float64,Int64},P)
  q[:] = convert(Vector{Float64},q)

  !isa(constraints, Array) && (constraints = [constraints])
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

function assemble!(model::QOCS.Model,P::Real,q::Real,constraints::Union{QOCS.Constraint,Array{QOCS.Constraint}})
  Pm = spzeros(1,1)
  qm = zeros(1)
  Pm[1,1] = convert(Float64,P)
  qm[1] = convert(Float64,q)
  assemble!(model,Pm,qm,constraints)
end

assemble!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::AbstractMatrix{<:Real},constraints::Union{QOCS.Constraint,Array{QOCS.Constraint}}) = assemble!(model,P,vec(q),constraints)
assemble!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::Real,constraints::Union{QOCS.Constraint,Array{QOCS.Constraint}}) = assemble!(model,P,[q],constraints)

"""
    warmStart!(model,[x0,y0])

Provides the `QOCS.Model` with warm starting values for the primal variable `x` and/or the dual variable `y`.
"""
function warmStart!(model::QOCS.Model; x0::Union{Vector{Float64}, Nothing} = nothing, y0::Union{Vector{Float64}, Nothing} = nothing)
    if x0 isa Vector{Float64}
      if size(model.A,2) == length(x0)
        model.x0 = x0
      else
        error("Dimension of x0 doesn't match the dimension of A.")
      end
    end
    if y0 isa Vector{Float64}
      if length(model.b) == length(y0)
        model.y0 = y0
      else
        error("Dimension of y0 doesn't match the dimensions of the constraint variables A,b.")
      end
    end
end


"""
    set!(model,P,q,A,b,convexSets)

Sets model data directly based on provided fields.
"""
function set!(model::QOCS.Model,P::AbstractMatrix{<:Real},q::AbstractVector{<:Real},A::AbstractMatrix{<:Real},b::AbstractVector{<:Real},convexSets::Array{QOCS.AbstractConvexSet})
  # convert inputs
  P[:,:] = convert(SparseMatrixCSC{Float64,Int64},P)
  A[:,:] = convert(SparseMatrixCSC{Float64,Int64},A)
  q[:] = convert(Vector{Float64},q)
  b[:] = convert(Vector{Float64},b)

  model.P = P
  model.q = q
  model.A = A
  model.b = b
  model.m = length(b)
  model.n = length(q)
  model.convexSets = convexSets
  nothing
end
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
