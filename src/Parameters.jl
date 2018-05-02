module Parameters
using OSSDPTypes, Residuals, KKT
export setRhoVec!, adaptRhoVec!, updateRhoVec!


# set initial values of rhoVec
  function setRhoVec!(p::OSSDPTypes.Problem,settings::OSSDPTypes.OSSDPSettings)
    nEQ = p.K.f
    nINEQ = p.m - p.K.f

    p.ρVec = [1e3*settings.rho*ones(nEQ);settings.rho*ones(nINEQ)]
    p.Info.rho_updates[1] = settings.rho
    return nothing
  end


  # adapt rhoVec based on residual ratio
  function adaptRhoVec!(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
    # compute normalized residuals based on the working variables (dont unscale)
    IGNORE_SCALING = true
    r_prim, r_dual = calculateResiduals(ws,settings,IGNORE_SCALING)
    maxNormPrim, maxNormDual = maxResComponentNorm(ws,settings,IGNORE_SCALING)
    r_prim = r_prim/(maxNormPrim + 1e-10)
    r_dual = r_dual/(maxNormDual + 1e-10)

    newRho = settings.rho * sqrt(r_prim/(r_dual+1e-10))
    newRho = minimum([maximum([newRho,settings.RHO_MIN]),settings.RHO_MAX])
    # only update rho if significantly different than current rho
    # FIXME: Should it be settings.rho or previous rho???
    if (newRho > settings.adaptive_rho_tolerance*settings.rho) || (newRho < (1./settings.adaptive_rho_tolerance)*settings.rho)
      updateRhoVec!(newRho,ws.p,settings)
    end
    return nothing
  end

  function updateRhoVec!(newRho::Float64,p::OSSDPTypes.Problem,settings::OSSDPTypes.OSSDPSettings)
    nEQ = p.K.f
    nINEQ = p.m - p.K.f

    settings.rho = newRho
    p.ρVec = [1e3*newRho*ones(nEQ);newRho*ones(nINEQ)]
    # log rho updates to info variable
    push!(p.Info.rho_updates,newRho)
    factorKKT!(p,settings)
    return nothing
  end

end #MODULE