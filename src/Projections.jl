module Projections

using ..QOCS, LinearAlgebra
export nonNegativeOrthant!, zeroCone!,  freeCone!, box!, secondOrderCone!, sdcone!, projectCompositeCone!

# -------------------------------------
# Standard Projection Functions
# -------------------------------------

    function projectCompositeCone!(x::Vector{Float64},convexSets::Array{AbstractConvexSet})

      for convexSet in convexSets
        xpart = view(x,convexSet.indices)
        convexSet.project!(xpart,convexSet)
      end
      return x
    end


    # projection onto nonegative orthant R_+^n
    function nonNegativeOrthant!(x::SubArray{Float64},convexSet::QOCS.Nonnegatives)
      for i in eachindex(x)
        x[i] = max(x[i], 0.)
      end
      nothing
    end

    # projection onto zero cone
    function zeroCone!(x::SubArray{Float64},convexSet::QOCS.Zeros)
      for i in eachindex(x)
        x[i] = 0.
      end
      nothing
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

function floorsqrt!(s::Array,floor::AbstractFloat)
    @.s  = sqrt(max(floor,s))
end


end #module
