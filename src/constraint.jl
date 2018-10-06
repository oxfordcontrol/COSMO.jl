
"""
Constraint(A,b,convexSet,dim=0,indices=0:0)

Creates a QOCS constraint: `Ax + b ∈ convexSet`.

By default the following convex sets are supported: `Zeros`, `Nonnegatives`, `SecondOrderCone`, `PositiveSemidefiniteCone`.

# Examples
```jldoctest
julia> Constraint([1 0;0 1],zeros(2),QOCS.PositiveSemidefiniteCone())
Constraint
Size of A: (2, 2)
ConvexSet: QOCS.PositiveSemidefiniteCone
```

---
The optinal arguments `dim` and `indices` can be used to specify A and b for subparts of variable `x`. If `x` has dimension `dim=4`,
then x[2] and x[3] can be constrained to the zero cone in the following way:


# Examples
```jldoctest
julia> c = Constraint([1 0;0 1],zeros(2),QOCS.Zeros(),4,2:3)
Constraint
Size of A: (2, 4)
ConvexSet: QOCS.Zeros
```
Notice that extra columns of A have been added automatically.
```
julia>Matrix(c.A)
2×4 Array{Float64,2}:
0.0  1.0  0.0  0.0
0.0  0.0  1.0  0.0
```
"""
struct Constraint
    A::Union{AbstractMatrix{<:Real},AbstractVector{<:Real}}
    b::AbstractVector{<:Real}
    convexSet::AbstractConvexSet

    function Constraint(A::AbstractMatrix{<:Real},b::AbstractVector{<:Real},convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0)
        size(A,1) != length(b) && error("The dimensions of matrix A and vector b don't match.")
        size(b,2) != 1 && error("Input b must be a vector or a scalar.")
        if dim != 0
            dim < 0 && error("The dimension of x has to be greater than zero.")
        end

        A = convert(Matrix{Float64},A)
        b = convert(Vector{Float64},b)

        if indices != 0:0
            (indices.start < 1 || indices.stop < indices.start) && error("The index range for x has to be increasing and nonnegative.")
            dim < indices.stop && error("The dimension of x: $(dim) must be equal or higher than the the stop value of indices: $(indices.stop).")
            Ac = spzeros(size(A,1),dim)
            bc = zeros(dim)
            Ac[:,indices] = A
            bc[indices] = b
            A = Ac
            b = bc
        end
        convexSet.dim = size(A,1)
        new(A,b,convexSet)
    end

    function Constraint(A::Real,b::Real,convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0)
        Aarr = zeros(1,1)
        bvec = Vector{Float64}(undef,1)
        Aarr[1] = A
        bvec[1] = b
        Constraint(Aarr,bvec,convexSet,dim,indices)
    end

    Constraint(A::AbstractMatrix{<:Real},b::AbstractMatrix{<:Real},convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0) = Constraint(A,vec(b),convexSet,dim,indices)
    Constraint(A::AbstractMatrix{<:Real},b::Real,convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0) = Constraint(A,[b],convexSet,dim,indices)

end

function Base.show(io::IO, obj::QOCS.Constraint)
    print(io,"Constraint\nSize of A: $(size(obj.A))\nConvexSet: $(typeof(obj.convexSet))")
end
