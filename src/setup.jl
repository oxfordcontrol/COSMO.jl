function setup!(ws::COSMO.Workspace,settings::COSMO.Settings)
    # scale problem data
    if settings.scaling != 0
        scaleRuiz!(ws,settings)
    end

    setRhoVec!(ws,settings)

    # factor the KKT condition matrix
    ws.p.flags.FACTOR_LHS && factorKKT!(ws,settings)
end
