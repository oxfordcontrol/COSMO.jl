module Projections

using ..QOCS, LinearAlgebra
export nonNegativeOrthant!, zeroCone!,  freeCone!, box!, secondOrderCone!, sdcone!, projectCompositeCone!

# -------------------------------------
# Standard Projection Functions
# -------------------------------------

    function projectCompositeCone!(x::Vector{Float64},convexSets::Array{AbstractConvexSet},xprevious,k,sign)

      for convexSet in convexSets
        xpart = view(x,convexSet.indices)
        if isa(convexSet, QOCS.PositiveSemidefiniteCone)
          convexSet.project!(xpart, convexSet,xprevious,k,sign)
        else
          convexSet.project!(xpart,convexSet)
        end
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
      l = convexSet.l
      u = convexSet.u
      x[:] = min.( max.(x,l), u)
      nothing
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
        tNew = (normX+t)/2
        xt = (normX+t)/(2*normX).*xt
        x[:] = [tNew;xt]
      end
      nothing
    end

  # compute projection of X=mat(x) onto the positive semidefinite cone
   function sdcone!(x::SubArray{Float64},convexSet::QOCS.PositiveSemidefiniteCone)

    n = Int(sqrt(length(x)))

    # handle 1D case
    if size(x,1) == 1
      x = max.(x,0)
    else
      # recreate original matrix from input vectors
      #Xs = Symmetric(reshape(x,n,n))
      X = reshape(x,n,n)
      X = 0.5*(X+X')
      # compute eigenvalue decomposition
      F = eigen(X)

      ind = findall(x-> x>0, F.values)
      Λ = Matrix(Diagonal(F.values))
      UsE = F.vectors[:,ind]*sqrt.(Λ[ind,ind])
      Xp = UsE*UsE'
      # different method
      # Λ = diagm(F[:values])
      # Q = F[:vectors]
      # # set negative eigenvalues to 0
      # Xp = Q*max.(Λ,0)*Q'
      x[:] = vec(Xp)
    end
    nothing
  end

    # compute projection of X=mat(x) onto the positive semidefinite cone
    function sdcone!(x::SubArray{Float64},convexSet::QOCS.PositiveSemidefiniteCone,xprevious,k,sign)

      n = Int(sqrt(length(x)))
  
      # handle 1D case
      if size(x,1) == 1
        x = max.(x,0)
      else
        # recreate original matrix from input vectors
        #Xs = Symmetric(reshape(x,n,n))
        X = reshape(x,n,n)
        X = 0.5*(X+X')
        # compute eigenvalue decomposition
        F = eigen(X)
  
        ind = findall(x-> x>0, F.values)
        Λ = Matrix(Diagonal(F.values))
        UsE = F.vectors[:,ind]*sqrt.(Λ[ind,ind])
        Xp = UsE*UsE'
        # different method
        # Λ = diagm(F[:values])
        # Q = F[:vectors]
        # # set negative eigenvalues to 0
        # Xp = Q*max.(Λ,0)*Q'
        x[:] = vec(Xp)
      end
      nothing
    end



end #module