function setup!(ws::COSMO.Workspace)
	# scale problem data
	if ws.settings.scaling != 0
		scale_ruiz!(ws)
	end

	set_rho_vec!(ws)

	# factor the KKT condition matrix
	ws.flags.FACTOR_LHS && factor_KKT!(ws)
end
