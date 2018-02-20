module Residuals
  using OSSDPTypes
  export calculateResiduals, maxResComponentNorm, hasConverged

  function calculateResiduals(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
        n = ws.p.n
        m = ws.p.m
        H = [ws.p.A spzeros(m,n); speye(n) -speye(n)]
        u = [ws.x;ws.s]
        if settings.scaling != 0
          EinvAug = [ws.sm.Einv spzeros(m,n); spzeros(n,m) ws.sm.D]
          r_prim = norm((EinvAug*(H*u-[ws.p.b;spzeros(n)]))./ws.sm.sb,2)
          # ∇f0 + ∑ νi ∇hi(x*) == 0 condition
          r_dual = norm((ws.sm.cinv*ws.sm.Dinv*(ws.p.P*ws.x + ws.p.q + ws.λ + ws.p.A'*ws.ν))./ws.sm.sq,2)
        else
          r_prim = norm(H*u-[ws.p.b;spzeros(n)],2)
          r_dual = norm(ws.p.P*ws.x + ws.p.q + ws.λ + ws.p.A'*ws.ν,2)
        end

    return r_prim,r_dual
  end

  function maxResComponentNorm(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
    m = ws.p.m
    n = ws.p.n
    H = [ws.p.A spzeros(m,n); speye(n) -speye(n)]
    u = [ws.x;ws.s]
    if settings.scaling != 0
      EinvAug = [ws.sm.Einv spzeros(m,n); spzeros(n,m) ws.sm.D]
      maxNormPrim = max.(norm(EinvAug*H*u./ws.sm.sb,2), norm(sm.Einv*ws.p.b./ws.sm.sb,2))
      maxNormDual = max.(norm(ws.sm.cinv*ws.sm.Dinv*ws.p.P*ws.x./ws.sm.sq,2), norm(ws.sm.cinv*ws.sm.Dinv*ws.p.q./ws.sm.sq,2), norm(ws.sm.cinv*ws.sm.Dinv*ws.λ./ws.sm.sq,2), norm(ws.sm.cinv*ws.sm.Dinv*ws.p.A'*ws.ν./ws.sm.sq,2) )
    else
      maxNormPrim = max.(norm(H*u,2), norm(ws.p.b,2))
      maxNormDual = max.(norm(ws.p.P*ws.x,2), norm(ws.p.q,2), norm(ws.λ,2) ,norm(ws.p.A'*ws.ν,2) )
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