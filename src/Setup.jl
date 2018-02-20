module Setup
  using OSSDPTypes, Scaling, KKT, Parameters
  export setup!

  function setup!(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
    # scale problem data
    if settings.scaling != 0
      (settings.scaleFunc == 1) && scaleSCS!(ws,settings,K)
      (settings.scaleFunc == 2) && scaleRuiz!(ws,settings,K)
    end

    setRhoVec!(ws.p,settings)

    # factor the KKT condition matrix
    factorKKT!(ws.p,settings)
  end




end #END MODULE