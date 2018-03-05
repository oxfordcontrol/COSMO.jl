module Setup
  using OSSDPTypes, Scaling, KKT, Parameters
  export setup!

  function setup!(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
    # scale problem data
    if settings.scaling != 0
      (settings.scaleFunc == 1) && scaleSCS!(ws,settings)
      (settings.scaleFunc == 2) && scaleRuiz!(ws,settings)
    end
    setRhoVec!(ws.p,settings)

    # factor the KKT condition matrix
    factorKKT!(ws.p,settings)
  end




end #END MODULE
