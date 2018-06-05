module Infeasibility
using OSSDPTypes
export isPrimalInfeasible, isDualInfeasible

  # sup_{z in K_tilde_b = {-K} x {b} } <z,δy> = { <y,b> ,if y in Ktilde_polar
  #                                                 +∞   ,else}
  function supportFunction(δy,ws,settings)
    K = ws.p.K
    # check if δy is in polar cone to K_tilde = dual cone of K
    inPolarCone = true
    b = K.f + 1

    # check that relevant part of δy is in dual of R+
    if K.l > 0
      e = b + K.l - 1
      δy_sub = δy[b:e]
      inPolarCone = (minimum(δy_sub) >= -settings.eps_prim_inf)
      b = e+1
    end

    # check that relevant parts of δy is in dual of Lorenz cones
    if length(K.q) > 0 && inPolarCone
        for iii = 1:length(K.q)
          e = b + K.q[iii] - 1
          δy_sub = δy[b:e]
          inPolarCone = (norm(δy_sub[2:end],2) - δy_sub[1] <= settings.eps_prim_inf  )
          !inPolarCone && break
          b = e + 1
        end
    end

    # check that relevant parts of δy is in dual of SDP cones
    b = K.f + K.l + sum(K.q) + 1
    if length(K.s) > 0 && inPolarCone
        for iii = 1:length(K.s)
          e = b + K.s[iii] - 1
          δy_sub = δy[b:e]
          cDim = Int(sqrt(K.s[iii]))
          # FIXME: Here you might get complex eigenvalues, which causes problems with minimum(). Does it make sense that you can get complex eigenvalues in this problem type?
          # Current Fix: Just consider the real part
          inPolarCone = ( minimum(real(eig(reshape(full(δy_sub),cDim,cDim))[1])) >= -settings.eps_prim_inf)

          !inPolarCone && break
          b = e + 1
      end
    end

    if inPolarCone
      return (ws.p.b'*δy)[1]
    else
      return Inf
    end
  end

  function inRecessionCone(A_δx,ws,settings)
    K = ws.p.K

    inRecessionCone = true
    b = 1

    # check that relevant parts of Aδx are in recession cone of {0}
    if K.f  > 0
        e = b + K.f - 1
        Aδx_sub = A_δx[b:e]
        inRecessionCone = (norm(Aδx_sub,Inf) <= settings.eps_dual_inf)
        b = e + 1
    end

    # check that relevant parts of Aδx are in recession cone of R-
    if K.l > 0 && inRecessionCone
      e = b + K.l - 1
      Aδx_sub = A_δx[b:e]
      inRecessionCone = (maximum(Aδx_sub) <= settings.eps_dual_inf)
      b = e +1
    end

    # check that relevant parts of Aδx are in recession cone of negative Lorenz cones
    # i.e. check that ||x|| <= -t
    if length(K.q) > 0 && inRecessionCone
        for iii = 1:length(K.q)
          e = b + K.q[iii] - 1
          Aδx_sub = A_δx[b:e]
          inRecessionCone = (norm(Aδx_sub[2:end],2) + Aδx_sub[1] <= settings.eps_dual_inf)
          !inRecessionCone && break
          b = e + 1
        end
    end


    # check that relevant parts of Aδx are in recession cone of negative semidefinite cones
    b = K.f + K.l + sum(K.q) + 1
    if length(K.s) > 0 && inRecessionCone
        for iii = 1:length(K.s)
          e = b + K.s[iii] - 1
          Aδx_sub = A_δx[b:e]
          cDim = Int(sqrt(K.s[iii]))
          inRecessionCone = ( maximum(real(eig(reshape(full(Aδx_sub),cDim,cDim))[1])) <= settings.eps_dual_inf)
          !inRecessionCone && break
          b = e + 1
      end
    end

    return inRecessionCone
  end



  function isPrimalInfeasible(δy,ws,settings)
    # calculate unscaled norm of δy
    if settings.scaling != 0
      norm_δy = norm(ws.sm.E*δy,Inf)
    else
      norm_δy = norm(δy,Inf)
    end

    # make sure norm is unequal to zero before continuing
    if norm_δy > settings.eps_prim_inf
      # test condition A'δy = 0
      A_δy = ws.p.A'*δy
      settings.scaling != 0 && (A_δy = ws.sm.Dinv*A_δy)

      if norm(A_δy,Inf)/norm_δy <= settings.eps_prim_inf
        # test condition S_K(δy) < 0
        sF = supportFunction(δy/norm_δy,ws,settings)
        if sF <= settings.eps_prim_inf
           return true
        end
      end
    end
    return false
  end


  function isDualInfeasible(δx,ws,settings)
    # calculate unscaled norm of δx
    if settings.scaling != 0
      norm_δx = norm(ws.sm.D*δx,Inf)
    else
      norm_δx = norm(δx,Inf)
    end

    if norm_δx > settings.eps_dual_inf
      # test condition <q,δx> < 0
      if (ws.p.q'*δx)[1]/(norm_δx*ws.sm.c) < -settings.eps_dual_inf
        # test condition Pδx == 0
        P_δx = ws.p.P*δx
        settings.scaling != 0 && (P_δx = ws.sm.Dinv*P_δx)
        if norm(P_δx,Inf)/(norm_δx*ws.sm.c) <= settings.eps_dual_inf

          # test condition Ax in Ktilde_b∞
          A_δx = ws.p.A*δx
          settings.scaling != 0 && (A_δx = ws.sm.Einv*A_δx)
            if inRecessionCone(A_δx/norm_δx,ws,settings)
              return true
            end
          end
      end
    end
    return false
  end

end #MODULE