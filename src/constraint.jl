
"""
	Constraint{T <: AbstractFloat}(A, b, convex_set_type, dim = 0, indices = 0:0)

Creates a COSMO constraint: `Ax + b ∈ convex_set`.

By default the following convex set types are supported: `ZeroSet`, `Nonnegatives`, `SecondOrderCone`, `PsdCone`, `PsdConeTriangle`.

# Examples
```jldoctest; setup = :(using COSMO)
julia> COSMO.Constraint([1 0;0 1], zeros(2), COSMO.Nonnegatives)
Constraint
Size of A: (2, 2)
ConvexSet: COSMO.Nonnegatives{Float64}
```

For convex sets that require their own data, it is possible to pass the pass the instantiated object directly rather than the type name.
# Examples
```jldoctest; setup = :(using COSMO)
julia> COSMO.Constraint([1 0;0 1], zeros(2), COSMO.Box([-1.;-1.],[1.;1.]))
Constraint
Size of A: (2, 2)
ConvexSet: COSMO.Box{Float64}
```


---
The optional arguments `dim` and `indices` can be used to specify A and b for subparts of variable `x`. If `x` has dimension `dim = 4`,
then x[2] and x[3] can be constrained to the zero cone in the following way:


# Examples
```jldoctest; setup = :(using COSMO)
julia> c = COSMO.Constraint([1 0;0 1], zeros(2), COSMO.ZeroSet, 4, 2:3)
Constraint
Size of A: (2, 4)
ConvexSet: COSMO.ZeroSet{Float64}
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
		convex_set::AbstractConvexSet{T},
		dim::Integer = 0,
		indices::UnitRange = 0:0) where {T <: AbstractFloat}

		size(A, 1) != length(b) && throw(DimensionMismatch("The dimensions of matrix A and vector b don't match."))
		size(A,1)  != convex_set.dim && throw(DimensionMismatch("The row dimension of A doesn't match the dimension of the constraint set."))
		size(b, 2) != 1 && throw(DimensionMismatch("Input b must be a vector or a scalar."))

		if indices != 0:0
			(indices.start < 1 || indices.stop < indices.start) && throw(DomainError("The index range for x has to be increasing and nonnegative."))
			dim < indices.stop && throw(DomainError("The dimension of x: $(dim) must be equal or higher than the the stop value of indices: $(indices.stop)."))
			Ac = spzeros(T, size(A, 1), dim)
			Ac[:, indices] = A
			A = Ac
		end
		new(A, b, convex_set)
	end
end

#We always choose to internally convert whatever datatype
#we are given to a floating point type either specified by
#the user or to DefaultFloat
function Constraint(A::AbstractMatrix{T}, b::AbstractVector{T}, args...) where {T <: Real}
	if T <: AbstractFloat
		return Constraint{T}(A, b, args...)
	else
		return Constraint{DefaultFloat}(Base.convert(AbstractMatrix{DefaultFloat}, A), Base.convert(AbstractVector{DefaultFloat}, b), args...)
	end
 end

#support the case where only a Type is specified for the set, and we need to
#create an instance of the appropriate size
function Constraint{T}(
	A::AbstractMatrix{T},
	b::AbstractVector{T},
	set_type::Type{ <: AbstractConvexSet}, args...) where {T <: AbstractFloat}

	# this constructor doesnt work with cones that need special arguments like the power cone
	set_type <: ArgumentCones && error("You can't create a constraint by passing the convex set as a type, if your convex set is a $(set_type). Please pass an object.")
	# call the appropriate set constructor
    n = size(A, 1)
    # we can deduce whether the PsdCone must be real or complex from the dimension
	if set_type <: PsdConeTriangles
	    if n == 1 || isqrt(n)^2 != n
            convex_set = set_type{T, T}(n)
        else
            convex_set = set_type{T, Complex{T}}(n)
        end
    else
        convex_set = set_type{T}(n)
    end
	Constraint{T}(A, b, convex_set, args...)
end


#allow various ways of passing A and b and convert  to AbstractMatrix and AbstractVector types

#all others convert first to matrix / vector args
function Constraint(A::T, b::T, args...) where {T <: Real}
	Constraint(reshape([A], 1, 1), [b], args...)
end

function Constraint(A::AbstractMatrix{T}, b::AbstractMatrix{T}, args...) where {T <: Real}
	Constraint(A, vec(b), args...)
end

function Constraint(A::Union{AbstractVector{T}, AbstractMatrix{T}}, b::T, args...) where {T <: Real}
	Constraint(reshape(A, 1, length(A)), [b], args...)
end

function Constraint(A::AbstractVector{T}, b::AbstractVector{T}, args...) where {T <: Real}
	Constraint(reshape(A, length(A), 1), b, args...)
end


function Base.show(io::IO, obj::COSMO.Constraint{T}) where {T <: AbstractFloat}
	print(io,"Constraint{$(T)}\nSize of A: $(size(obj.A))\nConvexSet: $(typeof(obj.convex_set))")
end
