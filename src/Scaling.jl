module Scaling

using ..QOCS, SparseArrays, LinearAlgebra, Statistics
export scaleRuiz!,reverseScaling!


  function kktColNorms!(P,A,normLHS,normRHS)

    colNorms!(normLHS,P,reset = true);   #start from zero
    colNorms!(normLHS,A,reset = false);  #incrementally from P norms
    rowNorms!(normRHS,A)                 #same as column norms of A'
    return nothing
  end

  function limitScaling!(s::Array,set::QOCS.Settings)
      @.s = clip(s,set.MIN_SCALING,set.MAX_SCALING,1.)
      return nothing
  end
  function limitScaling(s::Number,set::QOCS.Settings)
      s = clip(s,set.MIN_SCALING,set.MAX_SCALING,1.)
      return s
  end



  function scaleRuiz!(ws::QOCS.Workspace,set::QOCS.Settings)

      #references to scaling matrices from workspace
      D    = ws.sm.D
      E    = ws.sm.E
      c    = ws.sm.c[]

      #unit scaling to start
      D.diag .= 1.
      E.diag .= 1.
      c       = 1.

      #use the inverse scalings as intermediate
      #work vectors as well, since we don't
      #compute the inverse scaling until the
      #final step
      Dwork = ws.sm.Dinv
      Ework = ws.sm.Einv

      #references to QP data matrices
      P = ws.p.P
      A = ws.p.A
      q = ws.p.q
      b = ws.p.b

      #perform scaling operations for a fixed
      #number of steps, i.e. no tolerance or
      #convergence check
      for i = 1:set.scaling

          kktColNorms!(P,A,Dwork.diag,Ework.diag)
          limitScaling!(Dwork.diag,set)
          limitScaling!(Ework.diag,set)

          invsqrt!(Dwork.diag)
          invsqrt!(Ework.diag)

          # Scale the problem data and update the
          # equilibration matrices
          scaleData!(P,A,q,b,Dwork,Ework,1.)
          lmul!(Dwork,D)        #D[:,:] = Dtemp*D
          lmul!(Ework,E)        #D[:,:] = Dtemp*D

          # now use the Dwork array to hold the
          # column norms of the newly scaled P
          # so that we can compute the mean
          colNorms!(Dwork.diag,P)
          mean_col_norm_P = mean(Dwork.diag)
          inf_norm_q      = norm(q,Inf)

          if mean_col_norm_P  != 0. && inf_norm_q != 0.

            inf_norm_q = limitScaling(inf_norm_q,set)
            scale_cost = max(inf_norm_q,mean_col_norm_P)
            scale_cost = limitScaling(scale_cost,set)
            ctmp = 1.0 / scale_cost

            # scale the penalty terms and overall scaling
            P[:]  *= ctmp
            q    .*= ctmp
            c     *= ctmp
          end

      end #end Ruiz scaling loop

      # for certain cones we can only use a
      # a single scalar value.  In these cases
      # compute an adjustment to the overall scaling
      # so that the aggregate scaling on the cone
      # in questions turns out to be component-wise eq
      if rectifySetScalings!(E,Ework,ws.p.convexSets)
          #only rescale if the above returns true,
          #i.e. some cone scalings were rectified
          scaleData!(P,A,q,b,I,Ework,1.)
-         lmul!(Ework,E)
     end

      #scale set components
      scaleSets!(E,ws.p.convexSets)

      #update the inverse scaling data, c and c_inv
      ws.sm.Dinv.diag .= 1. ./ D.diag
      ws.sm.Einv.diag .= 1. ./ E.diag

      #These are Base.RefValue type so that
      #scaling can remain an immutable
      ws.sm.c[]        = c
      ws.sm.cinv[]     = 1. ./ c

      # scale the potentially warm started variables
      ws.x[:] = ws.sm.Dinv * ws.x
      ws.μ[:] = ws.sm.Einv * ws.μ
      ws.μ  .*= c

      return nothing
  end

  function invsqrt!(A::Array{T}) where{T}
      @. A = oneunit(T) / sqrt(A)
  end

  function rectifySetScalings!(E,Ework,sets)

      anyRectifiedBlocks  = false
      Ework.diag         .= 1.

      # NB : we should actually provide each cone
      # with the opportunity to provide its own
      # (possibly non-scalar) rectification

      for set in sets
          isScalar, = set.scale!(set)

          @views if isScalar
              #at least one block was scalar
              anyRectifiedBlocks = true
              tmp = mean(E.diag[set.indices])
              Ework.diag[set.indices] .= tmp./E.diag[set.indices]
          end
      end
      return anyRectifiedBlocks
  end


  function scaleSets!(E,sets)

      # scale set components (like u,l in a box)
      for set in sets
          scaleInfo = set.scale!(set)
          if length(scaleInfo) > 1
              #NB : Memory allocated here?
              for elem in scaleInfo[2:end]
                  elem[set.indices] .*= E.diag[set.indices]
              end
          end
      end
  end


  function scaleData!(P,A,q,b,Ds,Es,cs=1.)

      lrmul!(Ds,P,Ds) # P[:,:] = Ds*P*Ds
      lrmul!(Es,A,Ds) # A[:,:] = Es*A*Ds
      q[:] = Ds*q
      b[:] = Es*b

      if cs != 1.
          P .*= cs
          q .*= cs
      end
      return nothing
  end


  function reverseScaling!(ws::QOCS.Workspace)

    cinv = ws.sm.cinv[] #from the Base.RefValue type

    ws.x[:] = ws.sm.D*ws.x
    ws.s[:] = ws.sm.Einv*ws.s
    ws.ν[:] = ws.sm.E*ws.ν
    ws.μ[:] = ws.sm.E*ws.μ
    ws.ν  .*= cinv
    ws.μ  .*= cinv

    # reverse scaling for model data
    if ws.p.flags.REVERSE_SCALE_PROBLEM_DATA
        scaleData!(ws.p.P,ws.p.A,ws.p.q,ws.p.b,
                   ws.sm.Dinv,ws.sm.Einv,cinv)
    end
    return nothing
  end

end # MODULE
