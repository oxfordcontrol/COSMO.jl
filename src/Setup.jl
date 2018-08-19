module Setup
  using ..QOCS,..Scaling, ..KKT, ..Parameters
  export setup!

  function setup!(ws::QOCS.WorkSpace,settings::QOCS.Settings)
    # scale problem data
    if settings.scaling != 0
      scaleRuiz!(ws,settings)
    end

    setRhoVec!(ws,settings)

    # factor the KKT condition matrix
    ws.p.flags.FACTOR_LHS && factorKKT!(ws,settings)
  end




end #END MODULE
