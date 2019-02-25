function setup!(ws::COSMO.Workspace)

  allocate_set_memory!(ws)
	# scale problem data
	if ws.settings.scaling != 0
		scale_ruiz!(ws)
	end

	set_rho_vec!(ws)

	# factor the KKT condition matrix
	ws.flags.FACTOR_LHS && factor_KKT!(ws)
end

function allocate_set_memory!(ws)
  for (k, C) in enumerate(ws.p.C.sets)
    allocate_memory!(C)
  end
end