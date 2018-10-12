module Projections

using ..QOCS, LinearAlgebra
export project!, projectCompositeCone!

# -------------------------------------
# Standard Projection Functions
# -------------------------------------

    function projectCompositeCone!(s_views::Vector{SubArray},convexSets::Array{AbstractConvexSet})
      for (i,convexSet) in enumerate(convexSets)
        project!(s_views[i],convexSet)
      end
    end


    # projection onto nonegative orthant R_+^n
    function project!(x::SubArray{T},convexSet::QOCS.Nonnegatives) where{T}
      @.x = max(x,zero(T))
    end

    # projection onto zero cone
    function project!(x::SubArray{T},convexSet::QOCS.Zeros) where{T}
      x .= zero(T)
    end


    # compute projection of x onto a box defined by l and u
    function project!(x::SubArray{Float64},convexSet::QOCS.Box)
      l  = convexSet.l
      u  = convexSet.u
      x .= clip.(x,l,u)
    end


    # projection onto second-order-cone {(t,x) | ||x||_2 <= t}
    function project!(x::SubArray{Float64},convexSet::QOCS.SecondOrderCone)
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
   function project!(x::SubArray{Float64},convexSet::QOCS.PositiveSemidefiniteCone)

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

function floorsqrt!(s::Array,floor::AbstractFloat)
    @.s  = sqrt(max(floor,s))
end


end #module
