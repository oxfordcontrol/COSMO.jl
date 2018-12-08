function setup!(ws::COSMO.Workspace, settings::COSMO.Settings)
	# scale problem data
	if settings.scaling != 0
		scale_ruiz!(ws, settings)
	end

	set_rho_vec!(ws, settings)

	# factor the KKT condition matrix
	ws.p.flags.FACTOR_LHS && factor_KKT!(ws, settings)
end
