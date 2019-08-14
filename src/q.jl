mutable struct UpdatableQ{T} <: Factorization{T}
    """
    Gives the Q1 matrix of a qr factorization of an (n, m) matrix
    where Q := [Q1 Q2] and Q is an (n, m) orthornormal matrix
    """
    Q::Matrix{T}
    n::Int
    m::Int
    reset_counter::Int

    Q1::SubArray{T, 2, Matrix{T}, Tuple{Base.Slice{Base.OneTo{Int}}, UnitRange{Int}}, true}

    function UpdatableQ(A::AbstractMatrix{T}) where {T}
        n, m = size(A)
        @assert(m <= n, "Too many columns in the matrix.")

        F = qr(A)
        Q = zeros(T, n, n)
        Q[:, 1:m] .= Matrix(F.Q)

        new{T}(Q, n, m, 0, view(Q, :, 1:m))
    end
end

function set_Q!(F::UpdatableQ, Q::Matrix{T}) where {T}
  @assert(F.n == size(Q, 1))
  F.m = size(Q, 2)
  update_views!(F)
  copyto!(F.Q1, Q)
end

function add_column!(F::UpdatableQ{T}, a::AbstractVector{T}) where {T}
    @inbounds begin
        q_new = view(F.Q, :, F.m + 1)
        copyto!(q_new, a)
    end

    @inbounds for i in 1:F.m
        q = view(F.Q, :, i)
        r = dot(q, q_new)
        axpy!(-r, q, q_new)
    end
   
    d = norm(q_new)
    @. q_new ./= d

    if abs(d) > 1e-6
      F.m += 1; update_views!(F)
    else
        nothing
      # @warn "Skipping singular columng"
    end

    nothing
end

function update_views!(F::UpdatableQ{T}) where {T}
    F.Q1 = view(F.Q, :, 1:F.m)
end

function add_columns!(F::UpdatableQ{T}, A::AbstractMatrix{T}) where {T}
    if F.reset_counter < 10
        F.reset_counter += 1
        for i = 1:size(A, 2)
            add_column!(F, view(A, :, i))
        end
        return false
    else
        F.reset_counter = 0
        Q1 = Matrix(qr([F.Q1 A]).Q)
        # @show norm(F.Q1 - Q1[:, 1:F.m])
        F.m += size(A, 2)
        update_views!(F)
        F.Q1 .= Q1
        return true
    end
end
