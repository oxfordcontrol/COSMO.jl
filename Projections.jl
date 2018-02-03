module Projections

export nonNegativeOrthant, zeroCone,  freeCone, box, secondOrderCone, sdcone

# -------------------------------------
# HELPER FUNCTIONS
# -------------------------------------

    # projection onto nonegative orthant R_+^n
    function nonNegativeOrthant(x)
      return max.(x,0)
    end

    # projection onto zero cone
    function zeroCone(x)
      return zeros(size(x,1),1)
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
    function secondOrderCone(x,t::Float64)
      normX = norm(x,2)
      if  normX <= -t
        return 0.0.*x,0
      elseif normX <= t
        return x,t
      else
        tNew = (normX+t)/2
        x = (normX+t)/(2*normX).*x
        return x,tNew
      end
    end

  # compute projection of X=mat(x) onto the positive semidefinite cone
   function sdcone(x,n::Int64)
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