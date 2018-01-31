module Scaling

export scaleProblem!

  function normKKTCols(P,A)
    normPCols = [norm(P[:,i],Inf) for i in 1:size(P,2)]
    normACols = [norm(A[:,i],Inf) for i in 1:size(A,2)]
    normLeftPart = max.(normPCols,normACols)
    normATCols = [norm(A'[:,i],Inf) for i in 1:size(A,1)]

    return [normLeftPart;normATCols]
  end

  function limitScaling!(δVec,set)

    # Array case
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

    return δVec
  end

  function scaleProblem!(problem,scaleMatrices,set)
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
      δVec = limitScaling!(δVec,set)
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
      P = Dtemp*(P*Dtemp)
      A = Etemp*A*Dtemp
      q = Dtemp*q
      b = Etemp*b

      # Update equilibrium matrices D and E
      D = Dtemp*D
      E = Etemp*E

      # # Second step cost normalization
      # norm_P_cols = spla.norm(P, np.inf, axis=0).mean()
      # inf_norm_q = np.linalg.norm(q, np.inf)
      # inf_norm_q = self._limit_scaling(inf_norm_q)
      # scale_cost = np.maximum(inf_norm_q, norm_P_cols)
      # scale_cost = self._limit_scaling(scale_cost)
      # scale_cost = 1. / scale_cost
    end

    # write scaled problem data to object
    problem.P = P
    problem.A = A
    problem.b = b
    problem.q = q

    return D,E
  end


end # MODULE