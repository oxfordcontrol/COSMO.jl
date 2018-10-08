module Projections

using ..QOCS, LinearAlgebra
export nonNegativeOrthant!, zeroCone!,  freeCone!, box!, secondOrderCone!, sdcone!, projectCompositeCone!

# -------------------------------------
# Standard Projection Functions
# -------------------------------------

    function projectCompositeCone!(x::Vector{Float64},convexSets::Array{AbstractConvexSet}, use_lanczos::Bool)

      for convexSet in convexSets
        xpart = view(x,convexSet.indices)
        if use_lanczos && isa(convexSet, QOCS.PositiveSemidefiniteCone)
          sdcone_lanczos!(xpart,convexSet)
        else
          convexSet.project!(xpart,convexSet)
        end
      end
      return x
    end


    # projection onto nonegative orthant R_+^n
    function nonNegativeOrthant!(x::SubArray{T},convexSet::QOCS.Nonnegatives) where{T}
      @.x = max(x,zero(T))
    end

    # projection onto zero cone
    function zeroCone!(x::SubArray{T},convexSet::QOCS.Zeros) where{T}
      x .= zero(T)
    end


    # compute projection of x onto a box defined by l and u
    function box!(x::SubArray{Float64},convexSet::QOCS.Box)
      l  = convexSet.l
      u  = convexSet.u
      x .= clip.(x,l,u)
    end


    # projection onto second-order-cone {(t,x) | ||x||_2 <= t}
    function secondOrderCone!(x::SubArray{Float64},convexSet::QOCS.SecondOrderCone)
      t = x[1]
      xt = view(x,2:length(x))

      normX = norm(xt,2)
      if normX <= -t
        x[:] .= 0.
      elseif normX <= t
        nothing
      else
        x[1] = (normX+t)/2
        #x(2:end) assigned via view
        @. xt   = (normX+t)/(2*normX)*xt
      end
      nothing
    end

  # compute projection of X=mat(x) onto the positive semidefinite cone
   function sdcone!(x::SubArray{Float64},convexSet::QOCS.PositiveSemidefiniteCone)
    n = Int(sqrt(length(x)))

    # handle 1D case
    if size(x,1) == 1
      x = max.(x,0.)
    else
      # symmetrized square view of x
      X    = reshape(x,n,n)
      X[:] = 0.5*(X+X')
      # compute eigenvalue decomposition
      # then round eigs up and rebuild
      s,U  = eigen!(X)
      floorsqrt!(s,0.)
      rmul!(U,Diagonal(s))
      mul!(X, U, U')
    end
    nothing
  end

  function sdcone_lanczos!(x::SubArray{Float64},convexSet::QOCS.PositiveSemidefiniteCone)
    n = Int(sqrt(length(x)))
    m = convexSet.subspace_dimension
    X = reshape(x,n,n)
    X = (X+X')/2

    if m < 0 || m > n/2  # should be m >= n/k where k is number of lanczos steps
      # First iteration
      E = eigen(X)
      if sum(E.values .> 0.0) >= sum(E.values .< 0.0)
        convexSet.positive_subspace = true
      else
        convexSet.positive_subspace = false
      end
      ind = findall(x-> x>0, E.values)
      Z = E.vectors[:, ind]*Diagonal(E.values[ind])
      Xp = Z*E.vectors[:,ind]'
    else
      z = view(convexSet.subspace, 1:m*n)
      Z = reshape(z, n, m)
      if !convexSet.positive_subspace
        X .= -X;
      end
      # Just two steps of block-Lanczos
      Q = Array(qr([Z X*Z]).Q)
      QXQ = Q'*X*Q; QXQ = (QXQ+QXQ')/2;
      E = eigen(QXQ);
      ind = findall(x-> x>0, E.values)
      Z = Q*E.vectors[:, max(minimum(ind.-2),1):maximum(ind)]

      Xp = Q*E.vectors*max.(Diagonal(E.values),0)*E.vectors'*Q';
      Xp = (Xp+Xp')/2;

      if !convexSet.positive_subspace
          Xp .= -X .+ Xp;
      end
    end
    x[:] = vec(Xp)
    convexSet.subspace_dimension = size(Z)[2]
    convexSet.subspace[1:convexSet.subspace_dimension*n] = vec(Z)

    nothing
  end

function floorsqrt!(s::Array,floor::AbstractFloat)
    @.s  = sqrt(max(floor,s))
end


end #module
