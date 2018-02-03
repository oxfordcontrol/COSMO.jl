module Scaling

using OSSDPTypes
export scaleProblem!,reverseScaling!

  function normKKTCols(P::SparseMatrixCSC{Float64,Int64},A::SparseMatrixCSC{Float64,Int64})
    normPCols = [norm(P[:,i],Inf) for i in 1:size(P,2)]
    normACols = [norm(A[:,i],Inf) for i in 1:size(A,2)]
    normLeftPart = max.(normPCols,normACols)
    normATCols = [norm(A[i,:],Inf) for i in 1:size(A,1)]

    return [normLeftPart;normATCols]
  end

  function limitScaling!(δVec::Union{Float64,Array{Float64,1}},set::OSSDPTypes.sdpSettings)

    if length(δVec) > 1
      for iii = 1:length(δVec)
        if δVec[iii] < set.MIN_SCALING
          δVec[iii] = 1.0
        elseif δVec[iii] > set.MAX_SCALING
          δVec[iii] = set.MAX_SCALING
        end
      end
    # scalar case
    else
      if δVec < set.MIN_SCALING
        δVec = 1.0
      elseif δVec > set.MAX_SCALING
        δVec = set.MAX_SCALING
      end
    end
    return nothing
  end

  function scaleProblem!(problem::OSSDPTypes.problem,scaleMatrices::OSSDPTypes.scaleMatrices,set::OSSDPTypes.sdpSettings)
    P = problem.P
    A = problem.A
    b = problem.b
    q = problem.q
    m = problem.m
    n = problem.n

    c = 1.0
    sTemp = ones(n+m)

    #initialize scaling matrices
    D = speye(n)
    Dtemp = speye(n)
    E = speye(m)
    Etemp = speye(m)
    if m == 0
      E = 0
      Etemp = 0
    end


    for iii = 1:set.scaling

      # First step Ruiz
      δVec = normKKTCols(P,A)
      limitScaling!(δVec,set)
      δVec = sqrt(δVec)
      sTemp = 1./δVec

      # Obtain scaling matrices
      Dtemp = spdiagm(sTemp[1:n])
      if m == 0
        Etemp = 0
      else
        Etemp = spdiagm(sTemp[n+1:end])
      end

      # Scale data
      P[:,:] = Dtemp*(P*Dtemp)
      A[:,:] = Etemp*A*Dtemp
      q[:] = Dtemp*q
      b[:] = Etemp*b

      # Update equilibrium matrices D and E
      D = Dtemp*D
      E = Etemp*E

      # Second step cost normalization
      norm_P_cols = mean([norm(P[:,i],Inf) for i in 1:size(P,2)])
      inf_norm_q = norm(q,Inf)
      limitScaling!(inf_norm_q,set)
      scale_cost = maximum([inf_norm_q norm_P_cols])
      limitScaling!(scale_cost,set)
      scale_cost = 1. / scale_cost
      c_temp = scale_cost

      # Normalize cost
      P[:,:] = c_temp * P
      q[:] = c_temp * q

      # Update scaling
      c = c_temp * c


    end
    scaleMatrices.D = D
    scaleMatrices.E = E
    scaleMatrices.Dinv = spdiagm(1./diag(D))
    scaleMatrices.Einv = spdiagm(1./diag(E))
    scaleMatrices.c = c
    scaleMatrices.cinv = 1./c
    return D,E
  end


  function reverseScaling!(x::Array{Float64,1},s::Array{Float64,1},P::SparseMatrixCSC{Float64,Int64},q::SparseVector{Float64,Int64},sm::OSSDPTypes.scaleMatrices)
    x[:] = sm.D*x
    P[:,:] = sm.Dinv*P*sm.Dinv
    q[:] = sm.Dinv*q
    s[:] = sm.D*s
  end

end # MODULE