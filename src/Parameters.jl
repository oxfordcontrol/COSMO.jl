module Parameters
using ..QOCS,..Residuals, ..KKT
export setRhoVec!, adaptRhoVec!, updateRhoVec!


# set initial values of rhoVec
  function setRhoVec!(ws::QOCS.Workspace,settings::QOCS.Settings)
    p = ws.p
    # nEQ = p.K.f
    # nINEQ = p.m - p.K.f
    ws.ρ = settings.rho
    ws.ρVec = ws.ρ*ones(p.m) #[1e3*ws.ρ*ones(nEQ);ws.ρ*ones(nINEQ)]
    push!(ws.Info.rho_updates,ws.ρ)
    return nothing
  end


  # adapt rhoVec based on residual ratio
  function adaptRhoVec!(ws::QOCS.Workspace,settings::QOCS.Settings)
    # compute normalized residuals based on the working variables (dont unscale)
    IGNORE_SCALING = true
    r_prim::Float64, r_dual::Float64 = calculateResiduals(ws,settings,IGNORE_SCALING)
    maxNormPrim::Float64, maxNormDual::Float64  = maxResComponentNorm(ws,settings,IGNORE_SCALING)

    r_prim = r_prim/(maxNormPrim + 1e-10)
    r_dual = r_dual/(maxNormDual + 1e-10)

    newRho = ws.ρ * sqrt(r_prim/(r_dual+1e-10))
    newRho = min(max(newRho,settings.RHO_MIN),settings.RHO_MAX)
    # only update rho if significantly different than current rho
    # FIXME: Should it be settings.rho or previous rho???
    if (newRho > settings.adaptive_rho_tolerance*ws.ρ) || (newRho < (1 ./ settings.adaptive_rho_tolerance)*ws.ρ)
      updateRhoVec!(newRho,ws,settings)
    end
    return nothing
  end

  function updateRhoVec!(newRho::Float64,ws::QOCS.Workspace,settings::QOCS.Settings)
    p = ws.p
    nEQ = p.K.f
    nINEQ = p.m - p.K.f

    ws.ρ = newRho
    ws.ρVec = newRho*ones(p.m)#[1e3*newRho*ones(nEQ);newRho*ones(nINEQ)]
    # log rho updates to info variable
    push!(ws.Info.rho_updates,newRho)
    factorKKT!(ws,settings)
    return nothing
  end

end #MODULE
