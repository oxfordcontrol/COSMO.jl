module Projections

using ..QOCS, LinearAlgebra
export nonNegativeOrthant!, zeroCone!,  freeCone!, box!, secondOrderCone!, sdcone!, projectCompositeCone!
using Arpack

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
      V = E.vectors; λ = E.values;
      if sum(E.values .> 0.0) <= sum(E.values .< 0.0)
        convexSet.positive_subspace = true
        subspace_ind = λ .> 0
      else
        convexSet.positive_subspace = false
        subspace_ind = λ .< 0
      end
      Z = V[:, subspace_ind]
      ind = λ .> 0
      Xp = V[:, ind]*(λ[ind].*V[:, ind]')
    else
      z = view(convexSet.subspace, 1:m*n)
      Z = reshape(z, n, m)
      if !convexSet.positive_subspace
        X .= -X;
      end
      # Just one step of block-Lanczos
      Q = Array(qr([Z X*Z]).Q)
      QXQ = Q'*X*Q; QXQ = (QXQ+QXQ')/2;
      E = eigen(QXQ);
      l = E.values; V = E.vectors;
      # Sort eigenvalues-vectors
      sorted_idx = sortperm(l)
      l .= l[sorted_idx]; V .= V[:, sorted_idx];

      # Subspace to be used by next iterations
      idx = findfirst(l .> 0)
      Z = Q*V[:, max(idx-2,1):end]

      # Rayleigh values and vectors
      l_ = view(l, idx:length(l)); V_ = Q*view(V, :, idx:length(l))

      # Calculate projection error
      Xp = V_*(l_.*V_');
      if !convexSet.positive_subspace
          Xp .= -X .+ Xp;
      end

      # Calculate projection error
      P = I - V_*V_' # ToDo: Don't form the matrix
      λ_rem_max = eigs(Symmetric(P*X*P'), nev=1, ritzvec=false, which=:LR, tol=1e-4)[1]
      # λ_rem_max = sort(eigvals(Symmetric(P*X*P')))[end]

      projection_error = sqrt(
        2*norm(X*V_ - V_.*l_', 2)^2 +
        (n-size(V_, 2))*max(λ_rem_max[1], 0)^2
      )
      @show projection_error
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
