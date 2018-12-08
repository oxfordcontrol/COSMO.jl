
"""
Constraint(A, b, convex_set, dim = 0, indices = 0:0)

Creates a COSMO constraint: `Ax + b ∈ convex_set`.

By default the following convex sets are supported: `ZeroSet`, `Nonnegatives`, `SecondOrderCone`, `PsdCone`.

# Examples
```jldoctest
julia> Constraint([1 0;0 1], zeros(2), COSMO.PsdCone)
Constraint
Size of A: (2, 2)
ConvexSet: COSMO.PsdCone
```

---
The optinal arguments `dim` and `indices` can be used to specify A and b for subparts of variable `x`. If `x` has dimension `dim = 4`,
then x[2] and x[3] can be constrained to the zero cone in the following way:


# Examples
```jldoctest
julia> c = Constraint([1 0;0 1], zeros(2), COSMO.ZeroSet, 4, 2:3)
Constraint
Size of A: (2, 4)
ConvexSet: COSMO.ZeroSet
```
Notice that extra columns of A have been added automatically.
```
julia>Matrix(c.A)
2×4 Array{Float64,2}:
0.0  1.0  0.0  0.0
0.0  0.0  1.0  0.0
```
"""
struct Constraint{T <: AbstractFloat}
	A::Union{AbstractMatrix{T}, AbstractVector{T}}
	b::AbstractVector{T}
	convex_set::AbstractConvexSet{T}

	# constructor
	function Constraint{T}(
		A::AbstractMatrix{T},
		b::AbstractVector{T},
		set_type::Type{ <: AbstractConvexSet},
		dim::Integer = 0,
		indices::UnitRange = 0:0) where{T}

		size(A, 1) != length(b) && throw(DimensionMismatch("The dimensions of matrix A and vector b don't match."))
		size(b, 2) != 1 && throw(DimensionMismatch("Input b must be a vector or a scalar."))

		if indices != 0:0
			(indices.start < 1 || indices.stop < indices.start) && throw(DomainError("The index range for x has to be increasing and nonnegative."))
			dim < indices.stop && throw(DomainError("The dimension of x: $(dim) must be equal or higher than the the stop value of indices: $(indices.stop)."))
			Ac = spzeros(size(A, 1), dim)
			bc = zeros(dim)
			Ac[:, indices] = A
			bc[indices] = b
			A = Ac
			b = bc
		end
		dim = size(A, 1)
		# call the appropriate set constructor
		convex_set = set_type{T}(dim)
		new(A, b, convex_set)
	end
end

#We always choose to internally convert whatever datatype
#we are given to a floating point type either specified by
#the user or to DefaultFloat
Constraint(args...) = Constraint{DefaultFloat}(args...)

function Constraint{T}(A::AbstractMatrix, b::AbstractVector,
	set_type,
	dim::Integer = 0,
	indices::UnitRange = 0:0) where{T}
	Constraint{T}(AbstractMatrix{T}(A), AbstractVector{T}(b), set_type, dim, indices)
end

#all others convert first o matrix / vector args
function Constraint{T}(A::Real,b::Real,
	set_type,
	dim::Integer = 0,
	indices::UnitRange = 0:0) where{T}
	Constraint{T}(reshape([A], 1, 1), [b], set_type, dim, indices)
end

function Constraint{T}(A::AbstractMatrix,
	b::AbstractMatrix,
	set_type,
	dim::Integer = 0,
	indices::UnitRange = 0:0) where {T}

	Constraint{T}(A, vec(b), set_type, dim, indices)
end

function Constraint{T}(A::Union{AbstractVector,AbstractMatrix},
	b::Real,
	set_type,
	dim::Integer = 0,
	indices::UnitRange = 0:0) where{T}
	Constraint{T}(reshape(A, 1, length(A)), [b], set_type, dim, indices)
end


function Base.show(io::IO, obj::COSMO.Constraint)
	print(io,"Constraint\nSize of A: $(size(obj.A))\nConvexSet: $(typeof(obj.convex_set))")
end
