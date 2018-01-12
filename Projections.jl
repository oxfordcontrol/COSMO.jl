module Projections

export box, sdcone

# -------------------------------------
# HELPER FUNCTIONS
# -------------------------------------
  # compute projection of x onto a box defined by l and u
  function box(x::Array{Float64},l::Array{Float64},u::Array{Float64})
    return min.( max.(x,l), u)
  end

  # compute projection of X=mat(x) onto the positive semidefinite cone
   function sdcone(x::Array{Float64},n::Int64)
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