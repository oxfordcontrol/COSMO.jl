module Projections
using OSSDPTypes
export nonNegativeOrthant, zeroCone,  freeCone, box, secondOrderCone, sdcone, projectCompositeCone!

# -------------------------------------
# HELPER FUNCTIONS
# -------------------------------------

    function projectCompositeCone!(x,K::OSSDPTypes.Cone)
      b = 1

      if K.f  > 0
        e = b + K.f - 1
        zeroCone!(x, b, e)
        b = e + 1
      end

      if K.l > 0
        e = b + K.l - 1
        nonNegativeOrthant!(x, b, e)
        b = e +1
      end

      if length(K.q) > 0
        for iii = 1:length(K.q)
          e = b + K.q[iii] - 1
          x[b:e] = secondOrderCone(x[b:e])
          b = e + 1
        end
      end

      if length(K.s) > 0
        for iii = 1:length(K.s)
          #FIXME: Make sure they work directly on the input data, no copying
          e = b + K.s[iii] - 1
          x[b:e] = sdcone(x[b:e])
          b = e + 1
        end
      end
      return x
    end


    # projection onto nonegative orthant R_+^n
    function nonNegativeOrthant!(x, b, e)
      for i=b:e
        x[i] = max(x[i], 0.)
      end
      nothing
    end

    # projection onto zero cone
    function zeroCone!(x, b, e)
      for i=b:e
        x[i] = 0.
      end
      nothing
    end

      # projection onto free cone R^n
    function freeCone(x)
      return x
    end

    # compute projection of x onto a box defined by l and u
    function box(x,l,u)
      return min.( max.(x,l), u)
    end


    # projection onto second-order-cone {(t,x) | ||x||_2 <= t}
    function secondOrderCone(x)
      # FIXME: make sure no weird by reference is applied here
      t = x[1]
      x = x[2:length(x)]

      normX = norm(x,2)
      if  normX <= -t
        return [0.;0.0.*x]
      elseif normX <= t
        return [t;x]
      else
        tNew = (normX+t)/2
        x = (normX+t)/(2*normX).*x
        return [tNew;x]
      end
    end

  # compute projection of X=mat(x) onto the positive semidefinite cone
   function sdcone(x)

    n = Int(sqrt(length(x)))

    # handle 1D case
    if size(x,1) == 1
      return max.(x,0)
    else
      # recreate original matrix from input vectors
      X = reshape(x,n,n)
      X = X./2
      X = X+X'

      # compute eigenvalue decomposition
      F = eigfact(X)

      ind = find(x-> x>0, F[:values])
      Λ = diagm(F[:values])
      UsE = F[:vectors][:,ind]*sqrt.(Λ[ind,ind])
      Xp = UsE*UsE'
      # different method
      # Λ = diagm(F[:values])
      # Q = F[:vectors]
      # # set negative eigenvalues to 0
      # Xp = Q*max.(Λ,0)*Q'
      return vec(Xp)
    end
  end



end #module