# set initial values of rhoVec
function setRhoVec!(ws::COSMO.Workspace,settings::COSMO.Settings)
    m = ws.p.model_size[1]
    # nEQ = p.K.f
    # nINEQ = p.m - p.K.f
    ws.ρ = settings.rho
    ws.ρVec = ws.ρ * ones(m)

    # scale ρ values that belong to equality constraints with a factor of 1e3
    set_ind = findall(x -> typeof(x) == COSMO.ZeroSet{Float64}, ws.p.C.sets)
    if length(set_ind) > 0
        row_ind = COSMO.get_set_indices(ws.p.C.sets)
        for (i, rows) in enumerate(row_ind[set_ind])
            ws.ρVec[rows] *= 1e3
        end
    end
    push!(ws.Info.rho_updates, ws.ρ)
    return nothing
end


# adapt rhoVec based on residual ratio
function adaptRhoVec!(ws::COSMO.Workspace,settings::COSMO.Settings)
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

function updateRhoVec!(newRho::Float64,ws::COSMO.Workspace,settings::COSMO.Settings)

    ws.ρ     = newRho
    ws.ρVec .= newRho

     # scale ρ values that belong to equality constraints with a factor of 1e3
    set_ind = findall(x -> typeof(x) == COSMO.ZeroSet{Float64}, ws.p.C.sets)
    if length(set_ind) > 0
        row_ind = COSMO.get_set_indices(ws.p.C.sets)
        for (i, rows) in enumerate(row_ind[set_ind])
            ws.ρVec[rows] *= 1e3
        end
    end

    # log rho updates to info variable
    push!(ws.Info.rho_updates, newRho)
    factorKKT!(ws, settings)
    return nothing
end
