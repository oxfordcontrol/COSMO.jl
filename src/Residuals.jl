module Residuals
  using OSSDPTypes
  export calculateResiduals, maxResComponentNorm, hasConverged

  function calculateResiduals(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
        if settings.scaling != 0
          r_prim = norm((ws.sm.Einv*(ws.p.A*ws.x+ws.s-ws.p.b))./ws.sm.sb,Inf)
          # ∇f0 + ∑ νi ∇hi(x*) == 0 condition
          r_dual = norm((ws.sm.cinv*ws.sm.Dinv*(ws.p.P*ws.x + ws.p.q -ws.p.A'*ws.μ))./ws.sm.sq,Inf)
        else
          r_prim = norm(ws.p.A*ws.x+ws.s-ws.p.b,Inf)
          r_dual = norm(ws.p.P*ws.x + ws.p.q -ws.p.A'*ws.μ,Inf)
        end
        # FIXME: Why is it -A'μ ?
    return r_prim,r_dual
  end

  function maxResComponentNorm(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
    if settings.scaling != 0
      maxNormPrim = max.(norm(ws.sm.Einv*ws.p.A*ws.x./ws.sm.sb,Inf),  norm(ws.sm.Einv*ws.s./ws.sm.sb,Inf), norm(ws.sm.Einv*ws.p.b./ws.sm.sb,Inf))
      maxNormDual = max.(norm(ws.sm.cinv*ws.sm.Dinv*ws.p.P*ws.x./ws.sm.sq,Inf), norm(ws.sm.cinv*ws.sm.Dinv*ws.p.q./ws.sm.sq,Inf), norm(ws.sm.cinv*ws.sm.Dinv*ws.p.A'*ws.μ./ws.sm.sq,Inf) )
    else
      maxNormPrim = max.(norm(ws.p.A*ws.x,Inf),norm(ws.s,2), norm(ws.p.b,Inf))
      maxNormDual = max.(norm(ws.p.P*ws.x,Inf), norm(ws.p.q,2),norm(ws.p.A'*ws.μ,Inf) )
    end
    return maxNormPrim, maxNormDual
  end

  function hasConverged(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings,r_prim::Float64,r_dual::Float64)
    maxNormPrim, maxNormDual = maxResComponentNorm(ws,settings)
    ϵ_prim = settings.eps_abs + settings.eps_rel * maxNormPrim
    ϵ_dual = settings.eps_abs + settings.eps_rel * maxNormDual
    return ( r_prim < ϵ_prim  && r_dual < ϵ_dual)
  end

end #MODULE