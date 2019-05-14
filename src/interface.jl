"""
	assemble!(model, P, q, constraint(s); [settings, x0, y0, s0])

Assembles a `COSMO.Model` with a cost function defind by `P` and `q`, and a number of `constraints`.

The positive semidefinite matrix `P` and vector `q` are used to specify the cost function of the optimization problem:

```
min   1/2 x'Px + q'x
s.t.  Ax + b ∈ C
```
`constraints` is a `COSMO.Constraint` or an array of `COSMO.Constraint` objects that are used to describe the constraints on `x`.

---
The optional keyword argument `settings` can be used to pass custom solver settings:

```julia
custom_settings = COSMO.Settings(verbose = true);
assemble!(model, P, q, constraints, settings = custom_settings)
```
---
The optional keyword arguments `x0`, `s0`, and `y0` can be used to provide the solver with warm starting values for the primal variable `x`, the primal slack variable `s` and the dual variable `y`.

```julia
x_0 = [1.0; 5.0; 3.0]
COSMO.assemble!(model, P, q, constraints, x0 = x_0)
```

"""
function assemble!(model::Model{T},
	P::AbstractMatrix,
	q::AbstractVector,
	constraints::Union{Constraint{T},Vector{Constraint{T}}}; settings::COSMO.Settings = COSMO.Settings(),
	x0::Union{Vector{T}, Nothing} = nothing, y0::Union{Vector{T}, Nothing} = nothing, s0::Union{Vector{T}, Nothing} = nothing) where{ T<: AbstractFloat}


	!isa(constraints, Array) && (constraints = [constraints])

	merge_constraints!(constraints)
	model.p.P = P
	model.p.q = q
	n = length(q)
	m = sum(map( x-> x.dim, map( x-> x.convex_set, constraints)))

	model.p.model_size = [m;n]

	model.p.A = spzeros(Float64, m, n)
	model.p.b = spzeros(Float64, m)

	check_dimensions(model.p.P, model.p.q, model.p.A, model.p.b)

	# merge and sort the constraint sets
	sort!(constraints, by = x-> sort_sets(x.convex_set))
	row_num = 1
	for con in constraints
		process_constraint!(model.p, row_num, con.A, con.b, con.convex_set, n)
		row_num += con.convex_set.dim
	end

	# save the convex sets inside the model as a composite set
	model.p.C = CompositeConvexSet(map( x-> x.convex_set, constraints))
	model.settings = settings
	model.vars = Variables{T}(m, n, model.p.C)

	# if user provided (full) warm starting variables, update model
	x0 != nothing && warm_start_primal!(model, x0)
	s0 != nothing && warm_start_slack!(model, s0)
	y0 != nothing && warm_start_dual!(model, y0)
	nothing
end

# function assemble!(model::COSMO.Model{T},
# 	P::Real,q::Real,
# 	constraints::Union{COSMO.Constraint{T},Array{COSMO.Constraint{T}}}; settings::COSMO.Settings = COSMO.Settings()) where {T}
# 	Pm = spzeros(T, 1, 1)
# 	qm = zeros(T, 1)
# 	Pm[1, 1] = convert(T, P)
# 	qm[1] = convert(T, q)
# 	assemble!(model, Pm, qm, constraints, settings = settings)
# end

# Handle case where q is a 2-dimensional array instead of a 1-dimensional array
assemble!(model::COSMO.Model{T}, P::AbstractMatrix{T}, q::AbstractMatrix{T}, args...; kwargs...) where {T} = assemble!(model, P, vec(q), args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::AbstractVector, q::AbstractArray, args...; kwargs...) where {T} = assemble!(model, P[:, :], q, args...; kwargs...)
# Handle 1-D cases
assemble!(model::COSMO.Model{T}, P::Real, q::Real, args...; kwargs...) where {T} = assemble!(model, [P], [q], args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::Real, q::Union{AbstractMatrix, AbstractVector}, args...; kwargs...) where {T} = assemble!(model, [P], q, args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::Union{AbstractMatrix, AbstractVector}, q::Real, args...; kwargs...) where {T} = assemble!(model, P, [q], args...; kwargs...)

"""
	empty_model!(model)

Resets all the fields of `model` to that of a model created with `COSMO.Model()` (apart from the settings).
"""
function empty_model!(model::COSMO.Model{T}) where {T}
	model.p = ProblemData{T}()
	model.sm = ScaleMatrices{T}()
	model.vars = Variables{T}(1, 1, model.p.C)
	model.ρ = zero(T)
	model.ρvec = T[]
	model.kkt_solver = nothing
	model.flags = Flags()
	model.Info = Info([zero(T)])
	model.times = ResultTimes()
	nothing
end


function _warm_start!(z::Vector{T}, z0::Vector{T}, ind::Union{UnitRange{Int64}, Nothing}) where {T}
		ind == nothing && (ind = 1:length(z))
		length(ind) != length(z0) && throw(DimensionMismatch("Dimension of warm starting vector doesn't match the length of index range ind."))
		z[ind] = z0
end

"""
	warm_start_primal!(model, x0, [ind])

Provides the `COSMO.Model` with warm starting values for the primal variable `x`. `ind` can be used to warm start certain components of `x`.
"""
warm_start_primal!(model::COSMO.Model{T}, x0::Vector{T}, ind::Union{UnitRange{Int64}, Nothing}) where {T} = _warm_start!(model.vars.x, x0, ind)
warm_start_primal!(model::COSMO.Model{T}, x0::Vector{T}) where {T} = warm_start_primal!(model, x0, nothing)
warm_start_primal!(model::COSMO.Model{T}, x0::Real, ind::Int64) where {T} = (model.vars.x[ind] = x0)


"""
	warm_start_slack!(model, s0, [ind])

Provides the `COSMO.Model` with warm starting values for the primal slack variable `s`. `ind` can be used to warm start certain components of `s`.
"""
warm_start_slack!(model::COSMO.Model{T}, s0::Vector{T}, ind::Union{UnitRange{Int64}, Nothing}) where {T} = _warm_start!(model.vars.s.data, s0, ind)
warm_start_slack!(model::COSMO.Model{T}, s0::Vector{T}) where {T} = warm_start_slack!(model, s0, nothing)
warm_start_slack!(model::COSMO.Model{T}, s0::Real, ind::Int64) where {T} = (model.vars.s.data[ind] = s0)

# Notice that the sign of the dual variable y is inverted here, since internally the dual variable μ = -y is used
"""
	warm_start_dual!(model, y0, [ind])

Provides the `COSMO.Model` with warm starting values for the dual variable `y`. `ind` can be used to warm start certain components of `y`.
"""
warm_start_dual!(model::COSMO.Model{T}, y0::Vector{T}, ind::Union{UnitRange{Int64}, Nothing}) where {T} = _warm_start!(model.vars.μ, -y0, ind)
warm_start_dual!(model::COSMO.Model{T}, y0::Vector{T}) where {T} = warm_start_dual!(model, y0, nothing)
warm_start_dual!(model::COSMO.Model{T}, y0::Real, ind::Int64) where {T} = (model.vars.μ[ind] = -y0)

"""
	set!(model, P, q, A, b, convex_sets, [settings])

Sets model data directly based on provided fields.
"""
function set!(model::COSMO.Model,
	P::AbstractMatrix{<:Real},
	q::AbstractVector{<:Real},
	A::AbstractMatrix{<:Real},
	b::AbstractVector{<:Real},
	convex_sets::Vector{<: COSMO.AbstractConvexSet{T}}, settings::COSMO.Settings = COSMO.Settings()) where{T}

	check_dimensions(P, q, A, b)

	# convert inputs and copy them
	P_c = convert_copy(P, SparseMatrixCSC{Float64, Int64})
	A_c = convert_copy(A, SparseMatrixCSC{Float64, Int64})
	q_c = convert_copy(q, Vector{Float64})
	b_c = convert_copy(b, Vector{Float64})


	n = length(q)
	m = length(b)
	model.p.P = P_c
	model.p.q = q_c
	model.p.A = A_c
	model.p.b = b_c
	model.p.model_size = [m; n]
	model.p.C = CompositeConvexSet(deepcopy(convex_sets))
	model.vars = Variables{T}(m, n, model.p.C)
 	model.settings = settings
	nothing
end

function check_dimensions(P, q, A, b)
	size(A, 1) != length(b) && throw(DimensionMismatch("The dimensions of matrix A and vector b don't match."))
	size(A, 2) != length(q) && throw(DimensionMismatch("The dimensions of matrix A and vector q don't match."))
	size(b, 2) != 1 && throw(DimensionMismatch("Input b must be a vector or a scalar."))
	size(P, 1) != length(q) && throw(DimensionMismatch("The dimensions of matrix P and vector q don't match."))
	nothing
end

function check_A_dim(A::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}, n::Int64)
	size(A, 2) != n && throw(DimensionMismatch("The dimensions of a matrix A (m x $(size(A, 2))) in one of the constraints is inconsistent with the dimension of P ($(n))."))
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
  (C <: PsdConeTriangle) && return 6
  return 6
end

# transform A*x + b in {0}, to A*x + s == b, s in {0}
function process_constraint!(p::COSMO.ProblemData, row_num::Int64, A::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}, b::AbstractVector{<:Real}, C::AbstractConvexSet, n::Int64)
	check_A_dim(A, n)
	s = row_num
	e = row_num + C.dim - 1
	p.A[s:e, :] = -A
	p.b[s:e, :] = b
end
