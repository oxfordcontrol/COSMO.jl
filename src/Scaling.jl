module Scaling

using ..QOCS, SparseArrays, LinearAlgebra, Statistics
import LinearAlgebra: lmul!, rmul!
export scaleRuiz!,reverseScaling!

  function colNorms!(normArr::Array{Float64,1},A::SparseMatrixCSC{Float64,Int64},N::Int64)
    for i=1:N
      normArr[i] = norm(view(A,:,i),Inf)
    end
  end

  function normKKTCols!(δVec::Array{Float64,1},P::SparseMatrixCSC{Float64,Int64},A::SparseMatrixCSC{Float64,Int64},AT::SparseMatrixCSC{Float64,Int64},normPCols::Array{Float64,1},normACols::Array{Float64,1},normATCols::Array{Float64,1},m::Int64,n::Int64)
    colNorms!(normPCols,P,n)
    colNorms!(normACols,A,n)
    normLeftPart = max.(normPCols,normACols)
    AT[:,:] = permutedims(A)
    colNorms!(normATCols,AT,m)

    δVec[:] = [normLeftPart;normATCols]
    nothing
  end



  function limitScaling!(δVec::Union{Float64,Array{Float64,1}},set::QOCS.Settings)
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

  function scaleRuiz!(ws::QOCS.Workspace,set::QOCS.Settings)
    P = copy(ws.p.P)
    A = copy(ws.p.A)
    q = copy(ws.p.q)
    m = ws.p.m
    n = ws.p.n
    AT = spzeros(n,m)
    c = 1.0
    sTemp = ones(n+m)
    δVec = zeros(m+n)
    convexSets = ws.p.convexSets

    #initialize scaling matrices
    D = sparse(1.0I,n,n)
    Dtemp = sparse(1.0I,n,n)
    E = sparse(1.0I,m,m)
    Etemp = sparse(1.0I,m,m)
    if m == 0
      E = 0
      Etemp = 0
    end
    normPCols = zeros(n)
    normPCols2 = zeros(n)
    normACols = zeros(n)
    normATCols = zeros(m)

    for iii = 1:set.scaling

      # First step Ruiz
      QOCS.Scaling.normKKTCols!(δVec,P,A,AT,normPCols,normACols,normATCols,m,n)
      QOCS.Scaling.limitScaling!(δVec,set)
      δVec = sqrt.(δVec)
      sTemp = 1.0 ./δVec

      # Obtain scaling matrices
      Dtemp = sparse(Diagonal(sTemp[1:n]))
      if m == 0
        Etemp = 0
      else
        Etemp = sparse(Diagonal(sTemp[n+1:n+m]))
      end

      # Scale data
      P[:,:] = Dtemp*(P*Dtemp)
      A[:,:] = Etemp*A*Dtemp
      q[:] = Dtemp*q

      # Update equilibrium matrices D and E
      D[:,:] = Dtemp*D
      E[:,:] = Etemp*E
      # Second step cost normalization
      #colNorms!(normPCols2,P,n)
      #norm_P_cols = mean(normPCols2)
      norm_P_cols = mean([norm(P[:,i],Inf) for i in 1:size(P,2)])
      inf_norm_q = norm(q,Inf)
      if norm_P_cols != 0. && inf_norm_q != 0.
        QOCS.Scaling.limitScaling!(inf_norm_q,set)
        scale_cost = maximum([inf_norm_q norm_P_cols])
        QOCS.Scaling.limitScaling!(scale_cost,set)
        scale_cost = 1.0 ./ scale_cost
        c_temp = scale_cost

        # Normalize cost
        P[:,:] = c_temp * P
        q[:] = c_temp * q

        # Update scaling
        c = c_temp * c
      end
    end #END-Ruiz-Loop

    # make sure cone membership is preserved
    sTemp = [diag(D);diag(E)]


    for set in convexSets
      isScaleScalar, = set.scale!(set)
      if isScaleScalar
        ind = set.indices .+n

        sTemp[ind] .= mean(sTemp[ind])
      end
    end

    # Obtain scaling matrices
    D = sparse(Diagonal(sTemp[1:n]))
    if m == 0
      E = 0
    else
      E = sparse(Diagonal(sTemp[n+1:end]))
    end

    # perform final scaling
    ws.p.P = c*D*(ws.p.P*D)
    ws.p.A = E*ws.p.A*D
    ws.p.q = c*D*ws.p.q
    ws.p.b = E*ws.p.b

    # scale set components (like u,l in a box)
    for set in convexSets
      scaleInfo = set.scale!(set)
      if length(scaleInfo) > 1
        scaleVars = scaleInfo[2:end]
        for elem in scaleVars
          elem[:] = E*elem
        end
      end
    end

    ws.sm.D = D
    ws.sm.E = E
    ws.sm.Dinv = sparse(Diagonal(1 ./diag(D)))
    ws.sm.Einv = sparse(Diagonal(1 ./diag(E)))
    ws.sm.c = c
    ws.sm.cinv = 1 ./c

    # scale the potentially warm started variables
    ws.x[:] = ws.sm.Dinv *ws.x
    ws.μ[:] = ws.sm.Einv*ws.μ *c

    return nothing
  end



  function reverseScaling!(ws::QOCS.Workspace)
    D = ws.sm.D
    E = ws.sm.E
    Dinv = ws.sm.Dinv
    Einv = ws.sm.Einv
    c = ws.sm.c

    ws.x[:] = D*ws.x
    ws.s[:] = Einv*ws.s

    ws.ν[:] = E*ws.ν ./c
    ws.μ[:] = E*ws.μ ./c

    # reverse scaling for model data
    if ws.p.flags.REVERSE_SCALE_PROBLEM_DATA
      ws.p.P[:,:] = Dinv*ws.p.P*Dinv ./c
      ws.p.q[:] = (Dinv*ws.p.q) ./c
      ws.p.A[:,:] = Einv*ws.p.A*Dinv
      ws.p.b[:,:] = Einv*ws.p.b
    end
    return nothing
  end

end # MODULE



function lmul!(L::Diagonal, M::SparseMatrixCSC)

    #NB : Same as:  @views M.nzval .*= D.diag[M.rowval]
    #but this way allocates no memory at all and
    #is marginally faster
    for i = 1:(M.colptr[end]-1)
        M.nzval[i] *= L.diag[M.rowval[i]]
    end
    M
end

function rmul!(M::SparseMatrixCSC,R::Diagonal)
    for i = 1:M.n
        for j = M.colptr[i]:(M.colptr[i+1]-1)
            M.nzval[j] *= R.diag[i]
        end
    end
    M
end

function lrmul!(L::Diagonal, M::SparseMatrixCSC, R::Diagonal)
    for i = 1:M.n
        for j = M.colptr[i]:(M.colptr[i+1]-1)
            M.nzval[j] *= L.diag[M.rowval[j]]*R.diag[i]
        end
    end
    M
end

function lrmul!(L::Diagonal, M::AbstractMatrix, R::Diagonal)
    lmul!(L,rmul!(M,R))
end
