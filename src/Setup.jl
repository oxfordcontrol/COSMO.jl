module Setup
  using QOCS,..Scaling, ..KKT, ..Parameters
  export setup!

  function setup!(ws::QOCS.WorkSpace,settings::QOCS.Settings)
    # scale problem data
    if settings.scaling != 0
      scaleRuiz!(ws,settings)

    end

    setRhoVec!(ws.p,settings)

    # factor the KKT condition matrix
    factorKKT!(ws.p,settings)
  end




end #END MODULE
