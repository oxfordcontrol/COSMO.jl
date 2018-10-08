function setup!(ws::QOCS.Workspace,settings::QOCS.Settings)
    # scale problem data
    if settings.scaling != 0
        scaleRuiz!(ws,settings)
    end

    setRhoVec!(ws,settings)

    # factor the KKT condition matrix
    ws.p.flags.FACTOR_LHS && factorKKT!(ws,settings)
end
