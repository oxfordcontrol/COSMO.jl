module Scaling

using OSSDPTypes
export scaleRuiz!,reverseScaling!, scaleSCS!, scaleSymmetric!, findCloseSymmetricScaling


  function normKKTCols(P::SparseMatrixCSC{Float64,Int64},A::SparseMatrixCSC{Float64,Int64})
      normPCols = [norm(P[:,i],Inf) for i in 1:size(P,2)]
    normACols = [norm(A[:,i],Inf) for i in 1:size(A,2)]
    normLeftPart = max.(normPCols,normACols)
    normATCols = [norm(A[i,:],Inf) for i in 1:size(A,1)]

    return [normLeftPart;normATCols]
  end

  function TwonormKKTCols(P::SparseMatrixCSC{Float64,Int64},A::SparseMatrixCSC{Float64,Int64})
      normPCols = [norm(P[:,i],2) for i in 1:size(P,2)]
    normACols = [norm(A[:,i],2) for i in 1:size(A,2)]
    normLeftPart = max.(normPCols,normACols)
    normATCols = [norm(A[i,:],2) for i in 1:size(A,1)]

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
      δVec = sqrt.(δVec)
      sTemp = 1./δVec

      ix = n

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
          sTemp[ix+1:ix+numConeElem] = set.avgFunc(sTemp[ix+1:ix+numConeElem])
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
      if norm_P_cols != 0. && inf_norm_q != 0.
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

    end
    ws.sm.D = D
    ws.sm.E = E
    ws.sm.Dinv = spdiagm(1./diag(D))
    ws.sm.Einv = spdiagm(1./diag(E))
    ws.sm.c = c
    ws.sm.cinv = 1./c
    return nothing
  end

    # find cone and symmetry preserving scaling D*A*D that is close to Sm*A
  function findCloseSymmetricScaling(Sm::Array{Float64})

    SSm = Sm+Sm'

    # find max eigenvector of SSm and corresponding eigenvalue
    EVC = eig(SSm)
    λMax, maxInd = findmax(EVC[1])
    vMax = EVC[2][:,maxInd]
    d = sqrt(λMax/2)*vMax


    return diag(kron(diagm(d),diagm(d)))


  end
  function scaleSymmetric!(ws::OSSDPTypes.WorkSpace,set::OSSDPTypes.OSSDPSettings)
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
      δVec = sqrt.(δVec)
      sTemp = 1./δVec
      # start index as n, since the first n elem of sTemp define D and the remaining m define E (which has to preserve the cone)
      ix = n

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

      # variables in sdp cones. Use symmetric scalind D kron D to preserve cone membership
      if length(K.s) > 0
        for iii = 1:length(K.s)
          numConeElem = K.s[iii]
          matrixDim = Int(sqrt(numConeElem))
          Sm = reshape(sTemp[ix+1:ix+numConeElem],matrixDim,matrixDim)

          dVec = findCloseSymmetricScaling(Sm)


          sTemp[ix+1:ix+numConeElem] = dVec

          # perform a bunch of tests here to check correct implementation, i.e. distance to ideal scaling, preserving the cone membership, write some smaller functions to make tests easier

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
      if norm_P_cols != 0. && inf_norm_q != 0.
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

    end
    ws.sm.D = D
    ws.sm.E = E
    ws.sm.Dinv = spdiagm(1./diag(D))
    ws.sm.Einv = spdiagm(1./diag(E))
    ws.sm.c = c
    ws.sm.cinv = 1./c
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

    # find optimal scaling matrix for rows of M = [P A'
    #                                              A 0]
    # δVec is n+m x 1 with first n being the elements on diag of D and n+1:n+m being the elements on diag of E
    sTemp = TwonormKKTCols(P,A)

    # since scaling with E changes the rows of A (and also of s, since E*Ax+ E*s = E*b) we have to take care to preserve cone membership
    ix = n

    K.f > 0 && (ix += K.f)
    K.l > 0 && (ix += K.l)

    # # variables in second order cones. Use averaging to preserve cone membership
    # if length(K.q) > 0
    #   for iii = 1:length(K.q)
    #     numConeElem = K.q[iii]
    #     sTemp[ix+1:ix+numConeElem] = set.avgFunc(sTemp[ix+1:ix+numConeElem])
    #     ix += numConeElem
    #   end
    # end

    # # variables in sdp cones. Use averaging to preserve cone membership
    # if length(K.s) > 0
    #   for iii = 1:length(K.s)
    #     numConeElem = K.s[iii]
    #     sTemp[ix+1:ix+numConeElem] = set.avgFunc(sTemp[ix+1:ix+numConeElem])
    #     ix += numConeElem
    #   end
    # end

    limitScaling!(sTemp,set)
    sTemp = 1./sTemp

    D = spdiagm(sTemp[1:n])
    if m == 0
      E = 0
    else
      E = spdiagm(sTemp[n+1:end])
    end

    # Scale data
    P[:,:] = D*(P*D)
    A[:,:] = E*A*D
    q[:] = D*q
    b[:] = E*b

    # # scale q
    # # FIXME: Might make sense to pick row norms of P and A or sth like that
    # normARows = [norm(A[i,:],2) for i in 1:size(A,1)]
    # normPCols = [norm(P[:,i],2) for i in 1:size(P,2)]
    # meanRowNorm_AP = mean(vcat(normARows,normPCols))

    # # according to SCS: try to get rows  of A and q to have similar norm (for this solver P has to be scaled then as well)
    # nq = norm(q,2)

    # c = meanRowNorm_AP / maximum([nq set.MIN_SCALING])
    # q[:] = c*q
    # P[:] = c*P
    c = 1

    ws.sm.D = D
    ws.sm.E = E
    ws.sm.Dinv = spdiagm(1./diag(D))
    ws.sm.Einv = spdiagm(1./diag(E))
    ws.sm.c = c
    ws.sm.cinv = 1./c
    return nothing
  end






  function reverseScaling!(ws::OSSDPTypes.WorkSpace)
    # TODO: Double check what really is necessary to scale back (save time)
    D = ws.sm.D
    E = ws.sm.E
    Dinv = ws.sm.Dinv
    Einv = ws.sm.Einv
    c = ws.sm.c

    ws.x[:] = D*ws.x
    ws.p.P[:,:] = Dinv*ws.p.P*Dinv./c
    ws.p.q[:] = (Dinv*ws.p.q)./c
    ws.s[:] = Einv*ws.s
    # FIXME: Double check what has to be multiplied by scaling factors and what not
    ws.ν[:] = E*ws.ν./c
    ws.μ[:] = E*ws.μ./c
    return nothing
  end

end # MODULE