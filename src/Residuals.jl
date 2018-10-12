module Residuals
  using ..QOCS, LinearAlgebra
  export calculateResiduals, maxResComponentNorm, hasConverged

  function calculateResiduals(ws::QOCS.Workspace,settings::QOCS.Settings, IGNORESCALING_FLAG::Bool=false)
        if (settings.scaling != 0 && !IGNORESCALING_FLAG)
          r_prim = norm(ws.sm.Einv*(ws.p.A*ws.x+ws.s-ws.p.b),Inf)
          # ∇f0 + ∑ νi ∇hi(x*) == 0 condition
          r_dual = norm(ws.sm.cinv[]*ws.sm.Dinv*(ws.p.P*ws.x + ws.p.q -ws.p.A'*ws.μ),Inf)
        end
        if (settings.scaling == 0 || IGNORESCALING_FLAG )
          r_prim = norm(ws.p.A*ws.x+ws.s-ws.p.b,Inf)
          r_dual = norm(ws.p.P*ws.x + ws.p.q -ws.p.A'*ws.μ,Inf)
        end
        # FIXME: Why is it -A'μ ?
    return r_prim,r_dual
  end

  function maxResComponentNorm(ws::QOCS.Workspace,settings::QOCS.Settings, IGNORESCALING_FLAG::Bool=false)
    if (settings.scaling != 0 && !IGNORESCALING_FLAG)
      maxNormPrim = max.(norm(ws.sm.Einv*(ws.p.A*ws.x),Inf),  norm(ws.sm.Einv*ws.s,Inf), norm(ws.sm.Einv*ws.p.b,Inf))
      maxNormDual = max.(norm(ws.sm.cinv[]*ws.sm.Dinv*(ws.p.P*ws.x),Inf),
        norm(ws.sm.cinv[]*(ws.sm.Dinv*ws.p.q),Inf), norm(ws.sm.cinv[]*(ws.sm.Dinv*(ws.μ'*ws.p.A)'),Inf) )
    end
    if (settings.scaling == 0 || IGNORESCALING_FLAG)
      maxNormPrim = max.(norm(ws.p.A*ws.x,Inf),norm(ws.s,Inf), norm(ws.p.b,Inf))
      maxNormDual = max.(norm(ws.p.P*ws.x,Inf), norm(ws.p.q,Inf),norm(ws.p.A'*ws.μ,Inf) )
    end
    return maxNormPrim, maxNormDual
  end

  function hasConverged(ws::QOCS.Workspace,settings::QOCS.Settings,r_prim::Float64,r_dual::Float64)
    maxNormPrim, maxNormDual = maxResComponentNorm(ws,settings)
    ϵ_prim = settings.eps_abs + settings.eps_rel * maxNormPrim
    ϵ_dual = settings.eps_abs + settings.eps_rel * maxNormDual

    # if an optimal objective value was specified for the problem check if current solution is within specified accuracy
    obj_trueFLAG = true
    if !isnan(settings.obj_true)
      currentCost = ws.sm.cinv[]*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]
      obj_trueFLAG = abs(settings.obj_true-currentCost) <= settings.obj_true_tol
    end

    return ( r_prim < ϵ_prim  && r_dual < ϵ_dual && obj_trueFLAG)
  end

end #MODULE
