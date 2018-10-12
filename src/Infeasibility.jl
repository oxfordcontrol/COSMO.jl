module Infeasibility
using ..QOCS, LinearAlgebra
export isPrimalInfeasible, isDualInfeasible

  # sup_{z in K_tilde_b = {-K} x {b} } <z,δy> = { <y,b> ,if y in Ktilde_polar
  #                                                 +∞   ,else}

  function inDual(x,convexSets,tol)
    for convexSet in convexSets
      xpart = view(x,convexSet.indices)
      if !QOCS.inDual(xpart,convexSet,tol)
        return false
      end
    end
    return true
  end

  function supportFunction(δy,ws,settings)

    if inDual(δy,ws.p.convexSets,settings.eps_prim_inf)
      return (ws.p.b'*δy)[1]
    else
      return Inf
    end

  end

  function inRecessionCone(A_δx,ws,settings)
     convexSets = ws.p.convexSets
     for convexSet in convexSets
        A_δxpart = view(A_δx,convexSet.indices)
        if !QOCS.inRecc(A_δxpart,convexSet,settings.eps_dual_inf)
          return false
        end
      end
      return true
  end


  function isPrimalInfeasible(δy,ws,settings::Settings)
    # calculate unscaled norm of δy
    norm_δy = scalednorm(ws.sm.E,δy,Inf)::eltype(δy)

    # make sure norm is unequal to zero before continuing
    if norm_δy > settings.eps_prim_inf

      # test condition A'δy = 0
      A_δy = ws.p.A'*δy
      A_δy = ws.sm.Dinv*A_δy

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
    norm_δx = scalednorm(ws.sm.D,δx,Inf)::eltype(δx)

    if norm_δx > settings.eps_dual_inf
      # test condition <q,δx> < 0
      if  dot(ws.p.q,δx) / (norm_δx*ws.sm.c[]) < -settings.eps_dual_inf
        # test condition Pδx == 0
        P_δx = ws.p.P*δx
        P_δx = ws.sm.Dinv*P_δx
        if norm(P_δx,Inf)/(norm_δx*ws.sm.c[]) <= settings.eps_dual_inf

          # test condition Ax in Ktilde_b∞
          A_δx = ws.p.A*δx
          A_δx = ws.sm.Einv*A_δx
            if inRecessionCone(A_δx/norm_δx,ws,settings)
              return true
            end
          end
      end
    end
    return false
  end

end #MODULE
