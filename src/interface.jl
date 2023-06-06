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
The optional keyword arguments `x0` and `y0` can be used to provide the solver with warm starting values for the primal variable `x` and the dual variable `y`.

```julia
x_0 = [1.0; 5.0; 3.0]
COSMO.assemble!(model, P, q, constraints, x0 = x_0)
```

"""
function assemble!(model::Model{T},
	P::AbstractMatrix{T},
	q::AbstractVector{T},
	constraints::Union{Constraint{T}, Vector{Constraint{T}}}; settings::COSMO.Settings = COSMO.Settings{T}(),
	x0::Union{Vector{T}, Nothing} = nothing, y0::Union{Vector{T}, Nothing} = nothing) where { T <: AbstractFloat}


	!isa(constraints, Array) && (constraints = [constraints])
	eltype(settings) == T || throw(ArgumentError("The precision types of the model and the settings don't match."))
	type_checks(constraints)

	merge_constraints!(constraints)
	model.p.P = issparse(P) ? deepcopy(P) : sparse(P)
	model.p.q = deepcopy(q)
	n = length(q)
	m = sum(map( x-> x.dim, map( x-> x.convex_set, constraints)))

	model.p.model_size = [m; n]

	model.p.A = spzeros(T, m, n)
	model.p.b = spzeros(T, m)

	check_dimensions(model.p.P, model.p.q, model.p.A, model.p.b)

	# merge and sort the constraint sets
	sort!(constraints, by = x-> sort_sets(x.convex_set))
	row_num = 1
	for con in constraints
		process_constraint!(model.p, row_num, con.A, con.b, con.convex_set, n)
		row_num += con.convex_set.dim
	end

	# save the convex sets inside the model as a composite set
	model.p.C = CompositeConvexSet{T}(map( x-> x.convex_set, constraints))
	model.settings = deepcopy(settings)

	# the size of the temporary variables might change if the problem is decomposed
	# only allocate if it's not a cd problem
	pre_allocate_variables!(model)

	# set state of the model
	model.states.IS_ASSEMBLED = true

	# if user provided (full) warm starting variables, update model
	x0 != nothing && warm_start_primal!(model, x0)
	y0 != nothing && warm_start_dual!(model, y0)
	nothing
end


# Handle case where q is a 2-dimensional array instead of a 1-dimensional array
assemble!(model::COSMO.Model{T}, P::AbstractMatrix, q::AbstractMatrix, args...; kwargs...) where {T <: AbstractFloat} = assemble!(model, P, vec(q), args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::AbstractVector, q::AbstractArray, args...; kwargs...) where {T <: AbstractFloat} = assemble!(model, P[:, :], q, args...; kwargs...)
# Handle 1-D cases
assemble!(model::COSMO.Model{T}, P::Real, q::Real, args...; kwargs...) where {T <: AbstractFloat} = assemble!(model, Base.convert(T, P), Base.convert(T, q), args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::T, q::T, args...; kwargs...) where {T <: AbstractFloat} = assemble!(model, [P], [q], args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::Real, q::Union{AbstractMatrix{<: Real}, AbstractVector{<: Real}}, args...; kwargs...) where {T <: AbstractFloat} = assemble!(model, [P], q, args...; kwargs...)
assemble!(model::COSMO.Model{T}, P::Union{AbstractMatrix{<: Real}, AbstractVector{<: Real}}, q::Real, args...; kwargs...) where {T <: AbstractFloat} = assemble!(model, P, [q], args...; kwargs...)

# convert P, q to correct type
assemble!(model::COSMO.Model{T}, P::AbstractMatrix{Tp}, q::AbstractVector{Tq}, args...; kwargs...) where {T <: AbstractFloat, Tp <: Real, Tq <: Real} = assemble!(model, Base.convert(AbstractMatrix{T}, P), Base.convert(AbstractVector{T}, q), args...; kwargs...)
# Make sure constraints and model types are consistent
assemble!(model::COSMO.Model{T}, P::AbstractMatrix{T}, q::AbstractVector{T}, constraints::Union{Constraint{Tc}, Vector{Constraint{Tc}}}; kwargs...) where {T <: AbstractFloat, Tc <: Real} = throw(ArgumentError("The precision types of the model and the costraint(s) don't match."))
"""
	empty_model!(model)

Resets all the fields of `model` to that of a model created with `COSMO.Model()` (apart from the settings).
"""
function empty_model!(model::COSMO.Model{T}) where {T <: AbstractFloat}
	model.p = ProblemData{T}()
	model.sm = ScaleMatrices{T}()
	model.ci = ChordalInfo{T}()
	model.vars = Variables{T}(1, 1, model.p.C)
	model.utility_vars = UtilityVariables{T}(1, 1)
	model.ρ = zero(T)
	model.ρvec = T[]
	model.kkt_solver = nothing
	model.accelerator = CA.EmptyAccelerator()
	model.states = States()
	model.rho_updates = T[]
	model.rho_update_due = false
	model.times = ResultTimes()
	model.row_ranges = [0:0]
	nothing
end


function _warm_start!(z::AbstractVector{T}, z0::AbstractVector{T}, ind::Union{UnitRange{Int}, Nothing}) where {T <: AbstractFloat}
		ind == nothing && (ind = 1:length(z))
		length(ind) != length(z0) && throw(DimensionMismatch("Dimension of warm starting vector doesn't match the length of index range ind."))
		z[ind] = z0
end

"""
	warm_start_primal!(model, x0, [ind])

Provides the `COSMO.Model` with warm starting values for the primal variable `x`. `ind` can be used to warm start certain components of `x`.
"""
warm_start_primal!(model::COSMO.Model{T}, x0::Vector{T}, ind::Union{UnitRange{Int}, Nothing}) where {T <: AbstractFloat} = _warm_start!(model.vars.x, x0, ind)

warm_start_primal!(model::COSMO.Model{T}, x0::T, ind::Integer) where {T} = (model.vars.x[ind] = x0)

# if the full vector for x is provided, we can automatically warm start s = b - A*x as well
function warm_start_primal!(model::COSMO.Model{T}, x0::AbstractVector{T}) where {T <: AbstractFloat}
	warm_start_primal!(model, x0, nothing)

	# as scaling of (x, s, y) happens automatically in the solver, compute an un-scaled version of s0
	if model.states.IS_SCALED
		# scaled x0s
		mul!(model.utility_vars.vec_n,  model.sm.Dinv, x0)
		# scaled s0 = b - A * x0s
		mul!(model.utility_vars.vec_m,model.p.A, model.utility_vars.vec_n) # A * x0s
		@. model.utility_vars.vec_m *= -one(T)
		@. model.utility_vars.vec_m += model.p.b
		# unscaled s0 = Einv * s0s
		s0 = model.sm.Einv * model.utility_vars.vec_m
	else
		s0 = model.p.b - model.p.A * x0
	end
	warm_start_slack!(model, s0)
end

"""
	warm_start_slack!(model, s0, [ind])

Provides the `COSMO.Model` with warm starting values for the primal slack variable `s`. `ind` can be used to warm start certain components of `s`.
"""
warm_start_slack!(model::COSMO.Model{T}, s0::Vector{T}, ind::Union{UnitRange{Int}, Nothing}) where {T <: AbstractFloat} = _warm_start!(model.vars.s.data, s0, ind)
warm_start_slack!(model::COSMO.Model{T}, s0::Vector{T}) where {T} = warm_start_slack!(model, s0, nothing)
warm_start_slack!(model::COSMO.Model{T}, s0::T, ind::Integer) where {T} = (model.vars.s.data[ind] = s0)

# Notice that the sign of the dual variable y is inverted here, since internally the dual variable μ = -y is used
"""
	warm_start_dual!(model, y0, [ind])

Provides the `COSMO.Model` with warm starting values for the dual variable `y`. `ind` can be used to warm start certain components of `y`.
"""
warm_start_dual!(model::COSMO.Model{T}, y0::Vector{T}, ind::Union{UnitRange{Int}, Nothing}) where {T <: AbstractFloat} = _warm_start!(model.vars.μ, -y0, ind)
warm_start_dual!(model::COSMO.Model{T}, y0::Vector{T}) where {T} = warm_start_dual!(model, y0, nothing)
warm_start_dual!(model::COSMO.Model{T}, y0::T, ind::Integer) where {T} = (model.vars.μ[ind] = -y0)

"""
	warm_start!(model, x0, y0)

Provides the `COSMO.Model` with warm starting values for the primal variable `x` and the dual variable `y`.
"""
function warm_start!(model::COSMO.Model{T}, x0::Vector{T}, y0::Vector{T}) where {T <: AbstractFloat}
	warm_start_primal!(model, x0)
	warm_start_dual!(model, y0)
end


"""
	update!(model, q = nothing, b = nothing)

Updates the model data for `b` or `q`. This can be done without refactoring the KKT matrix. The vectors will be appropriatly scaled.
"""
function update!(model::COSMO.Model{T}; q::Union{Vector{T}, Nothing} = nothing, b::Union{Vector{T}, Nothing} = nothing) where {T <: AbstractFloat}
	m, n = model.p.model_size
	!model.states.IS_ASSEMBLED && error("Model has to be assembled once before one can start updating q or b.")
	if q != nothing
		length(q) != n && throw(DimensionMismatch("The dimension of q, does not agree with the model dimension, n."))
		model.states.IS_CHORDAL_DECOMPOSED && error("Problem vector q can not be updated if the model has been chordally decomposed before.")
		if model.states.IS_SCALED
			# if the internal model is scaled, also scale q
			mul!(model.p.q, model.sm.D, q)
			model.p.q .*= model.sm.c[]
		else
			@. model.p.q = q
		end
	end

	if b != nothing
		length(b) != m && throw(DimensionMismatch("The dimension of b, does not agree with the model dimension, m."))
		model.states.IS_CHORDAL_DECOMPOSED && error("Problem vector b can not be updated if the model has been chordally decomposed before.")
		if model.states.IS_SCALED
			mul!(model.p.b, model.sm.E, b)
		else
			@. model.p.b = b
		end
	end
end

"""
	set!(model, P, q, A, b, convex_sets, [settings])

Sets model data directly based on provided fields.
"""
function set!(model::COSMO.Model{T},
	P::AbstractMatrix{T},
	q::AbstractVector{T},
	A::AbstractMatrix{T},
	b::AbstractVector{T},
	convex_sets::Vector{<: COSMO.AbstractConvexSet{T}}, settings::COSMO.Settings{T} = COSMO.Settings{T}()) where {T <: AbstractFloat}

	check_dimensions(P, q, A, b)
	type_checks(convex_sets)

	# convert inputs and copy them
	P_c = convert_copy(P, SparseMatrixCSC{T, Int})
	A_c = convert_copy(A, SparseMatrixCSC{T, Int})
	q_c = convert_copy(q, Vector{T})
	b_c = convert_copy(b, Vector{T})


	n = length(q)
	m = length(b)
	model.p.P = P_c
	model.p.q = q_c
	model.p.A = A_c
	model.p.b = b_c
	model.p.model_size = [m; n]
	model.p.C = CompositeConvexSet{T}(deepcopy(convex_sets))

	pre_allocate_variables!(model)
 	model.settings = deepcopy(settings)

	# set state of the model
	model.states.IS_ASSEMBLED = true
	nothing
end

# a specific function that takes the sparse matrices P, A in (rowval, colptr, nzval)-form for easy interoperability with python interface
function set!(model::COSMO.Model{Tf},
	Prowval::Vector{Ti},
	Pcolptr::Vector{Ti},
	Pnzval::Vector{Tf},
	q::Vector{Tf},
	Arowval::Vector{Ti},
	Acolptr::Vector{Ti},
	Anzval::Vector{Tf},
	b::Vector{Tf},
	cone::Dict, l::Union{Nothing, Vector{Tf}}, u::Union{Nothing, Vector{Tf}}, m::Int, n::Int, settings::COSMO.Settings{Tf} = COSMO.Settings{Tf}()) where {Tf <: AbstractFloat, Ti <: Integer}

	# construct the sparse matrices
	if Ti == Int32
		Prowval = juliafy_integers(Prowval)
		Arowval = juliafy_integers(Arowval)
		Pcolptr = juliafy_integers(Pcolptr)
		Acolptr = juliafy_integers(Acolptr)
	end

	P = SparseMatrixCSC{Tf, Int}(n, n, Pcolptr, Prowval, Pnzval)
	A = SparseMatrixCSC{Tf, Int}(m, n, Acolptr, Arowval, Anzval)

	check_dimensions(P, q, A, b)

	model.p.P = P
	model.p.q = q
	model.p.A = A
	model.p.b = b
	model.p.model_size = [m; n]
	convex_sets = convex_sets_from_dict(cone, l, u)
	model.p.C = CompositeConvexSet{Tf}(convex_sets)

	pre_allocate_variables!(model)
 	model.settings = settings
	# set state of the model
	model.states.IS_ASSEMBLED = true

	nothing
end


# handle the case where settings is a transformed python dictionary
function set!(model::COSMO.Model{Tf},
	Prowval::Vector{Ti},
	Pcolptr::Vector{Ti},
	Pnzval::Vector{Tf},
	q::Vector{Tf},
	Arowval::Vector{Ti},
	Acolptr::Vector{Ti},
	Anzval::Vector{Tf},
	b::Vector{Tf},
	cone::Dict, l::Union{Nothing, Vector{Tf}}, u::Union{Nothing, Vector{Tf}}, m::Int, n::Int, settings_dict::Dict) where {Tf <: AbstractFloat, Ti <: Integer}

	settings = COSMO.Settings(settings_dict)
	COSMO.set!(model, Prowval, Pcolptr, Pnzval, q, Arowval, Acolptr, Anzval, b, cone, l, u, m, n, settings)

end

function juliafy_integers(arr::Vector{Int32})
	# 1-based indexing
	@. arr += 1
	# convert to 64bit
	return Base.convert.(Int, arr)
end

# given the cone-dict in scs format create an array of COSMO.AbstractConvexSet(s)
function convex_sets_from_dict(cone::Dict, l::Union{AbstractVector, Nothing}, u::Union{AbstractVector, Nothing})
	convex_sets = Vector{COSMO.AbstractConvexSet{Float64}}(undef, 0)
	haskey(cone, "f") && push!(convex_sets, COSMO.ZeroSet(cone["f"]))
	haskey(cone, "l") && push!(convex_sets, COSMO.Nonnegatives(cone["l"]))

	# second-order cones
	if haskey(cone, "q")
		socp_dim = cone["q"]
		for dim in socp_dim
			push!(convex_sets, COSMO.SecondOrderCone(dim))
		end
	end
	# sdp triangle cones
	if haskey(cone, "s")
		sdp_dim = cone["s"]
		for dim in sdp_dim
			push!(convex_sets, COSMO.PsdConeTriangle(dim))
		end
	end
	# primal exponential cones
	if haskey(cone, "ep")
		for k = 1:cone["ep"]
			push!(convex_sets, COSMO.ExponentialCone())
		end
	end
	# dual exponential cones
	if haskey(cone, "ed")
		for k = 1:cone["ed"]
			push!(convex_sets, COSMO.DualExponentialCone())
		end
	end
	# power cones
	if haskey(cone, "p")
		pow_exponents = cone["p"]
		for exponent in pow_exponents
			if exponent >= 0
				push!(convex_sets, COSMO.PowerCone(exponent))
			else
				push!(convex_sets, COSMO.DualPowerCone(-1. * exponent))
			end
		end
	end
	# box constraints
	if haskey(cone, "b")
		push!(convex_sets, COSMO.Box(l, u))
	end
	return convex_sets
end


function check_dimensions(P, q, A, b)
	size(A, 1) != length(b) && throw(DimensionMismatch("The dimensions of matrix A and vector b don't match."))
	size(A, 2) != length(q) && throw(DimensionMismatch("The dimensions of matrix A and vector q don't match."))
	size(b, 2) != 1 && throw(DimensionMismatch("Input b must be a vector or a scalar."))
	size(P, 1) != length(q) && throw(DimensionMismatch("The dimensions of matrix P and vector q don't match."))
	nothing
end

"Check whether the model will contain any PSD constraints with unsupported Floating-point precision."
function type_checks(convex_sets::Vector{<: COSMO.AbstractConvexSet{T}}) where {T <: AbstractFloat}
	for set in convex_sets
		type_checks(set)
	end
	return nothing
end
function type_checks(constraints::Vector{COSMO.Constraint{T}}) where {T <: AbstractFloat}
	for constraint in constraints
		type_checks(constraint.convex_set)
	end
	return nothing
end
type_checks(::AbstractConvexSet) = nothing
type_checks(::PsdCone{BigFloat}) = throw(ArgumentError("COSMO currently does not support the combination of PSD constraints and BigFloat."))
type_checks(::PsdConeTriangle{BigFloat}) = throw(ArgumentError("COSMO currently does not support the combination of PSD constraints and BigFloat."))



function check_A_dim(A::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}, n::Int)
	size(A, 2) != n && throw(DimensionMismatch("The dimensions of a matrix A (m x $(size(A, 2))) in one of the constraints is inconsistent with the dimension of P ($(n))."))
end

# convert x into type (which creates a copy) or copy x if type coincides
function convert_copy(x::AbstractArray, argtype::Type)
	if typeof(x) == argtype
		x_c = copy(x)
	else
		x_c = Base.convert(argtype, x)
	end
	return x_c
end

# merge zeros sets and nonnegative sets
function merge_constraints!(constraints::Array{COSMO.Constraint{T}}) where {T <: AbstractFloat}
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

function merge_zeros(constraints::Array{COSMO.Constraint{T}}) where {T <: AbstractFloat}
	m = sum(x -> x.dim, map(x -> x.convex_set, constraints))
	n = size(constraints[1].A, 2)
	A = spzeros(T, m, n)
	b = zeros(T, m)
	s = 1
	e = 0
	for cons in constraints
		e = s + cons.convex_set.dim - 1
		A[s:e, :] = cons.A
		b[s:e, :] = cons.b
		s = e + 1
	end
	return M = COSMO.Constraint{T}(A, b, ZeroSet)
end

function merge_nonnegatives(constraints::Array{COSMO.Constraint{T}}) where {T <: AbstractFloat}
	m = sum(x -> x.dim, map(x -> x.convex_set, constraints))
	n = size(constraints[1].A, 2)
	A = spzeros(T, m, n)
	b = zeros(T, m)

	s = 1
	e = 0
	for cons in constraints
		e = s + cons.convex_set.dim - 1
		A[s:e, :] = cons.A
		b[s:e, :] = cons.b
		s = e + 1
	end

	return M = COSMO.Constraint{T}(A, b, Nonnegatives)
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
function process_constraint!(p::COSMO.ProblemData{T}, row_num::Int, A::Union{AbstractVector{T}, AbstractMatrix{T}}, b::AbstractVector{T}, C::AbstractConvexSet{T}, n::Int) where {T <: AbstractFloat}
	check_A_dim(A, n)
	s = row_num
	e = row_num + C.dim - 1
	p.A[s:e, :] = -A
	p.b[s:e, :] = b
end


function pre_allocate_variables!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}
  m, n = ws.p.model_size
  ws.vars = Variables{T}(m, n, ws.p.C)
  ws.utility_vars = UtilityVariables{T}(m, n)
end
