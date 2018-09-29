module Scaling

using ..QOCS, SparseArrays, LinearAlgebra, Statistics
export scaleRuiz!,reverseScaling!


  function kktColNorms!(P,A,normLHS,normRHS)

    colNorms!(normLHS,P,reset = true);   #start from zero
    colNorms!(normLHS,A,reset = false);  #incrementally from P norms
    rowNorms!(normRHS,A)                 #same as column norms of A'
    return nothing
  end

  @inline function limitScaling!(s::Number,minval::Number,maxval::Number)
      s = s < minval ? 1  : (s > maxval ? maxval : s)
  end

  function limitScaling!(s::Array,minval,maxval)
      s .= limitScaling!.(s,minval,maxval)
  end

  function limitScaling!(s,set::QOCS.Settings)
      limitScaling!(s,set.MIN_SCALING,set.MAX_SCALING)
  end


  function scaleRuiz!(ws::QOCS.Workspace,set::QOCS.Settings)

      #references to scaling matrices from workspace
      D    = ws.sm.D
      E    = ws.sm.E
      DInv = ws.sm.Dinv
      EInv = ws.sm.Dinv

      #unit scaling to start
      D.diag .= 1.
      E.diag .= 1.
      c       = 1.

      #use the inverse scalings as intermediate
      #work vectors as well, since we don't
      #compute the inverse scaling until the
      #final step
      Dwork = Dinv
      Ework = Einv

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

          @. Dwork.diag = 1 / sqrt(Dwork.diag)
          @. Ework.diag = 1 / sqrt(Ework.diag)

          # Scale the problem data and update the
          # equilibration matrices
          scaleData(P,A,q,b,Dwork,Ework,1.)

          # Update equilibrium matrices D and E
          lmul!(Dwork,D)        #D[:,:] = Dtemp*D
          lmul!(Ework,E)        #D[:,:] = Dtemp*D

          # now use the Dwork array to hold the
          # column norms of the newly scaled P
          # so that we can compute the mean
          colNorms!(Dwork.diag,P)
          mean_col_norm_P = mean(Dwork.diag)
          inf_norm_q      = norm(q,Inf)

          if norm_P_cols != 0. && inf_norm_q != 0.

            limitScaling!(inf_norm_q,set)
            scale_cost = max(inf_norm_q,mean_col_norm_P)
            limitScaling!(scale_cost,set)
            c_tmp = 1.0 / scale_cost

            # scale the penalty terms and overall scaling
            P .*= c_temp
            q .*= c_temp
            c .*= c_temp
          end

      end #end Ruiz scaling loop

      # for certain cones we can only use a
      # a single scalar value.  In these cases
      # compute an adjustment to the overall scaling
      # so that the aggregate scaling on the cone
      # in questions turns out to be component-wise eq
      # NB : we should actually provide each cone
      # with the opportunity to provide its row
      # (possibly non-scalar) rectification

      anyScalarBlocks = false

      for set in convexSets
          isScalar, = set.scale!(set)

          if isScalar
              #at least one block was scalar
              anyScalarBlocks = true

              #NB: memory being allocated for ind here?
              ind = set.indices .+ n
              tmp = mean(D.diag[ind])
              Dwork.diag[ind] .= tmp./D.diag[ind]
          end

      end

      #if any adjustments were made, tweak the
      #scaling appropriately
      if anyScalarBlocks
          scaleData(P,A,q,b,Dwork,Ework,1.)
          # Update equilibrium matrices D and E
          lmul!(Dwork,D)        #D[:,:] = Dtemp*D
          lmul!(Ework,E)        #D[:,:] = Dtemp*D
      end

      #scale set components
      # scale set components (like u,l in a box)
      for set in convexSets
        scaleInfo = set.scale!(set)
        if length(scaleInfo) > 1
          #NB : Memory allocated here?
          scaleVars = scaleInfo[2:end]
          for elem in scaleVars
            elem[:] = E*elem
          end
        end
      end

      #update the inverse scaling data, c and c_inv
      ws.sm.Dinv.diag .= 1. ./ D.diag
      ws.sm.Einv.diag .= 1. ./ E.diag
      ws.sm.c          = c
      ws.sm.cinv       = 1. ./ c

      # scale the potentially warm started variables
      ws.x[:] = ws.sm.Dinv *ws.x
      ws.μ[:] = ws.sm.Einv*ws.μ *c

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

    ws.x[:] = ws.sm.D*ws.x
    ws.s[:] = ws.sm.Einv*ws.s
    ws.ν[:] = ws.sm.E*ws.ν*ws.sm.cinv
    ws.μ[:] = ws.sm.E*ws.μ*ws.sm.cinv

    # reverse scaling for model data
    if ws.p.flags.REVERSE_SCALE_PROBLEM_DATA
        scaleData!(ws.p.P,ws.p.A,ws.p.q,ws.p.b,
                   ws.sm.Dinv,ws.sm.Einv,ws.sm.cinv)
    end
    return nothing
  end

end # MODULE
