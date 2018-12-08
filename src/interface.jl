"""
assemble!(model, P, q, constraint(s), [x0, y0])

Assembles a `COSMO.Model` with a cost function defind by `P` and `q`, and a number of `constraints`.

The positive semidefinite matrix `P` and vector `q` are used to specify the cost function of the optimization problem:

```
min   1/2 x'Px + q'x
s.t.  Ax + b ∈ C
```
`constraints` is a `COSMO.Constraint` or an array of `COSMO.Constraint` objects that are used to describe the constraints on `x`.

---
The optinal arguments `x0` and `y0` can be used to provide the solver with warm starting values for the primal variable `x` and the dual variable `y`.
"""
function assemble!(model::Model{T},
	P::AbstractMatrix{T},
	q::AbstractVector{T},
	constraints::Union{Constraint{T},Vector{Constraint{T}}},
	x0::Union{Vector{Float64}, Nothing} = nothing, y0::Union{Vector{Float64}, Nothing} = nothing) where{ T<: AbstractFloat}

	# convert inputs
	#FIX ME : It should not be necessary to force sparsity,
	#since maybe we would like a dense solver.  Sparse for
	#now until we get a dense LDL option
	P_c = convert_copy(P, SparseMatrixCSC{Float64, Int64})
	q_c = convert_copy(q, Vector{Float64})

	!isa(constraints, Array) && (constraints = [constraints])
	# model.Flags.INFEASIBILITY_CHECKS = checkConstraintFunctions(constraints)

	merge_constraints!(constraints)
	model.P = P_c
	model.q = q_c
	n = length(q)
	m = sum(map( x-> x.dim, map( x-> x.convex_set, constraints)))

	model.model_size = [m;n]

	model.A = spzeros(Float64, m, n)
	model.b = spzeros(Float64, m)

	model.x0 = zeros(Float64, n)
	model.y0 = zeros(Float64, m)

	# merge and sort the constraint sets
	sort!(constraints, by = x-> sort_sets(x.convex_set))
	row_num = 1
	for con in constraints
		process_constraint!(model, row_num, con.A, con.b, con.convex_set)
		row_num += con.convex_set.dim
	end

	# save the convex sets inside the model as a composite set
	model.C = CompositeConvexSet(map( x-> x.convex_set, constraints))

	# if user provided warm starting variables, update model
	warm_start!(model, x0 = x0, y0 = y0)

	nothing
end

function assemble!(model::COSMO.Model,
	P::Real,q::Real,
	constraints::Union{COSMO.Constraint{T},Array{COSMO.Constraint{T}}}) where{T}
	Pm = spzeros(T, 1, 1)
	qm = zeros(T, 1)
	Pm[1, 1] = convert(T, P)
	qm[1] = convert(T, q)
	assemble!(model, Pm, qm, constraints)
end

function assemble!(model::COSMO.Model,
	P::AbstractMatrix{T},
	q::AbstractMatrix{T},
	constraints::Union{COSMO.Constraint{T},Array{COSMO.Constraint{T}}}) where{T}

	assemble!(model, P, vec(q), constraints)
end

assemble!(model::COSMO.Model, P::AbstractMatrix{<:Real}, q::Real, constraints::Union{COSMO.Constraint, Array{COSMO.Constraint}}) = assemble!(model, P, [q], constraints)

"""
warm_start!(model, [x0, y0])

Provides the `COSMO.Model` with warm starting values for the primal variable `x` and/or the dual variable `y`.
"""
function warm_start!(model::COSMO.Model; x0::Union{Vector{Float64}, Nothing} = nothing, y0::Union{Vector{Float64}, Nothing} = nothing)
	if x0 isa Vector{Float64}
		if size(model.A, 2) == length(x0)
			model.x0 = x0
		else
			error("Dimension of x0 doesn't match the dimension of A.")
		end
	end
	if y0 isa Vector{Float64}
		if length(model.b) == length(y0)
			model.y0 = y0
		else
			error("Dimension of y0 doesn't match the dimensions of the constraint variables A, b.")
		end
	end
end


"""
set!(model, P, q, A, b, convex_sets)

Sets model data directly based on provided fields.
"""
function set!(model::COSMO.Model,
	P::AbstractMatrix{<:Real},
	q::AbstractVector{<:Real},
	A::AbstractMatrix{<:Real},
	b::AbstractVector{<:Real},
	convex_sets::Vector{COSMO.AbstractConvexSet{T}}) where{T}

	# convert inputs and copy them
	P_c = convert_copy(P, SparseMatrixCSC{Float64, Int64})
	A_c = convert_copy(A, SparseMatrixCSC{Float64, Int64})
	q_c = convert_copy(q, Vector{Float64})
	b_c = convert_copy(b, Vector{Float64})

	model.P = P_c
	model.q = q_c
	model.A = A_c
	model.b = b_c
	model.model_size = [length(b); length(q)]
	model.C = CompositeConvexSet(convex_sets)
	nothing
end

# convert x into type (which creates a copy) or copy x if type coincides
function convert_copy(x::AbstractArray, type::Type)
	if typeof(x) == type
		x_c = copy(x)
	else
		x_c = convert(type, x)
	end
	return x_c
end

# merge zeros sets and nonnegative sets
function merge_constraints!(constraints::Array{COSMO.Constraint{T}}) where{T}
	# handle zeros sets
	ind = findall(set->typeof(set) == ZeroSet{T}, map(x -> x.convex_set, constraints))
	if length(ind) > 1
		M = merge_zeros(constraints[ind])
		deleteat!(constraints, ind)
		push!(constraints, M)
	end

	# handle nonnegative sets
	ind = findall(set->typeof(set) == Nonnegatives{T},map(x->x.convex_set,constraints))
	if length(ind) > 1
		M = merge_nonnegatives(constraints[ind])
		deleteat!(constraints, ind)
		push!(constraints, M)
	end
	nothing
end

function merge_zeros(constraints::Array{COSMO.Constraint{T}}) where{T}
	m = sum(x -> x.dim, map(x -> x.convex_set, constraints))
	n = size(constraints[1].A, 2)
	A = spzeros(m, n)
	b = zeros(m)
	s = 1
	e = 0
	for cons in constraints
		e = s + cons.convex_set.dim - 1
		A[s:e, :] = cons.A
		b[s:e, :] = cons.b
		s = e + 1
	end
	return M = COSMO.Constraint(A, b, ZeroSet)
end

function merge_nonnegatives(constraints::Array{COSMO.Constraint{T}}) where{T}
	m = sum(x -> x.dim, map(x -> x.convex_set, constraints))
	n = size(constraints[1].A, 2)
	A = spzeros(m, n)
	b = zeros(m)

	s = 1
	e = 0
	for cons in constraints
		e = s + cons.convex_set.dim - 1
		A[s:e, :] = cons.A
		b[s:e, :] = cons.b
		s = e + 1
	end

	return M = COSMO.Constraint(A, b, Nonnegatives)
end


function sort_sets(C::AbstractConvexSet)
	C = typeof(C)
	(C <: ZeroSet) && return 1
	(C <: Nonnegatives) && return 2
	(C <: Box) && return 3
	(C <: SecondOrderCone) && return 4
	(C <: PsdCone) && return 5
	return 6
end

# transform A*x + b in {0}, to A*x + s == b, s in {0}
function process_constraint!(model::COSMO.Model, row_num::Int64, A::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}, b::AbstractVector{<:Real}, C::AbstractConvexSet)
	s = row_num
	e = row_num + C.dim - 1
	model.A[s:e, :] = -A
	model.b[s:e, :] = b
end
