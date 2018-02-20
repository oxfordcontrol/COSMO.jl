module Scaling

using OSSDPTypes
export scaleRuiz!,reverseScaling!, scaleSCS!


  function normKKTCols(P::SparseMatrixCSC{Float64,Int64},A::SparseMatrixCSC{Float64,Int64})
    normPCols = [norm(P[:,i],Inf) for i in 1:size(P,2)]
    normACols = [norm(A[:,i],Inf) for i in 1:size(A,2)]
    normLeftPart = max.(normPCols,normACols)
    normATCols = [norm(A[i,:],Inf) for i in 1:size(A,1)]

    return [normLeftPart;normATCols]
  end

  function limitScaling!(δVec::Union{Float64,Array{Float64,1}},set::OSSDPTypes.OSSDPSettings)
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

  function scaleRuiz!(ws::OSSDPTypes.WorkSpace,set::OSSDPTypes.OSSDPSettings)
    P = ws.p.P
    A = ws.p.A
    b = ws.p.b
    q = ws.p.q
    m = ws.p.m
    n = ws.p.n
    K = ws.p.K
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
      ix = 0

      K.f > 0 && (ix += K.f)
      K.l > 0 && (ix += K.l)

      # variables in second order cones. Use averaging to preserve cone membership
      if length(K.q) > 0
        for iii = 1:length(K.q)
          numConeElem = K.q[iii]
          sTemp[ix+1:ix+numConeElem] = set.avgFunc(sTemp[ix+1:ix+numConeElem])
          ix += numConeElem
        end
      end

      # variables in sdp cones. Use averaging to preserve cone membership
      if length(K.s) > 0
        for iii = 1:length(K.s)
          numConeElem = K.s[iii]
          sTemp[ix+1:ix+numConeElem] = set.avgFunc(sTemp[ix+1:ix+numConeElem  ])
          ix += numConeElem
        end
      end

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
      D[:,:] = Dtemp*D
      E[:,:] = Etemp*E

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

#@show(D)
    end
 #   println("final D")
  #  @show(D)
    # @show(A)
    # @show(D)
    # @show(E)
    # @show(c)
    ws.sm.D = D
    ws.sm.E = E
    ws.sm.Dinv = spdiagm(1./diag(D))
    ws.sm.Einv = spdiagm(1./diag(E))
    ws.sm.c = c
    ws.sm.cinv = 1./c
    #@show(scaleMatrices)
    return nothing
  end

function scaleSCS!(ws::OSSDPTypes.WorkSpace,set::OSSDPTypes.OSSDPSettings)
    P = ws.p.P
    A = ws.p.A
    b = ws.p.b
    q = ws.p.q
    m = ws.p.m
    n = ws.p.n
    K = ws.p.K
    # scale rows of A' to have unit norm ()
    D =  [norm(A'[i,:],2) for i in 1:size(A',1)]
    # define row index of vector x
    ix = 0

    K.f > 0 && (ix += K.f)
    K.l > 0 && (ix += K.l)

    # variables in second order cones. Use averaging to preserve cone membership
    if length(K.q) > 0
      for iii = 1:length(K.q)
        numConeElem = K.q[iii]
        D[ix+1:ix+numConeElem] = set.avgFunc(D[ix+1:ix+numConeElem])
        ix += numConeElem
      end
    end

    # variables in sdp cones. Use averaging to preserve cone membership
    if length(K.s) > 0
      for iii = 1:length(K.s)
        numConeElem = K.s[iii]
        D[ix+1:ix+numConeElem] = set.avgFunc(D[ix+1:ix+numConeElem  ])
        ix += numConeElem
      end
    end
    limitScaling!(D,set)
    # scale rows by D first and then determine scaling for E
    D = spdiagm(1./D)
    A[:,:] = A*D

    E = [norm(A'[:,i],2) for i in 1:size(A',2)]
    limitScaling!(E,set)

    E = spdiagm(1./E)

    # Scale data on LHS
    P[:,:] = D*(P*D)
    A[:,:] = E*A

    # scale b
    # find mean row and col norm for scaled A
    meanRowNormA = mean([norm(A[i,:],2) for i in 1:size(A,1)])
    meanColNormA = mean([norm(A[:,i],2) for i in 1:size(A,2)])

    b[:] = E*b
    nb = norm(b,2)
    scale_b = meanColNormA / maximum([nb set.MIN_SCALING])
    #@show(scale_b)
    b[:] = scale_b * b
    #@show(b)

    # scale q
    q[:] = D*q
    nq = norm(q,2)
    scale_q = meanRowNormA / maximum([nq set.MIN_SCALING])
    #@show(scale_q)
    q[:] = scale_q*q
    #@show(q)

    ws.sm.D = D
    ws.sm.E = E
    ws.sm.Dinv = spdiagm(1./diag(D))
    ws.sm.Einv = spdiagm(1./diag(E))
    ws.sm.sq = scale_q
    ws.sm.sb = scale_b
    return nothing
  end




  function reverseScaling!(ws::OSSDPTypes.WorkSpace)
    # TODO: Double check what really is necessary to scale back (save time)
    D = ws.sm.D
    E = ws.sm.E
    Dinv = ws.sm.Dinv
    Einv = ws.sm.Einv
    c = ws.sm.c
    sb = ws.sm.sb
    sq = ws.sm.sq

    ws.x[:] = (D*ws.x)./sb
    ws.p.P[:,:] = Dinv*ws.p.P*Dinv/c
    ws.p.q[:] = (Dinv*ws.p.q)./sq./c
    ws.s[:] = (D*ws.s)./sb
    ws.ν[:] = E*ws.ν./c
    ws.λ[:] = Dinv*ws.λ./c
    return nothing
  end

end # MODULE