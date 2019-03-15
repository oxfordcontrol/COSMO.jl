function setup!(ws::COSMO.Workspace)
	# scale problem data
	if ws.settings.scaling != 0
		scale_ruiz!(ws)
	end

	set_rho_vec!(ws)

	# create a KKT solver object populated with our data
	if(ws.flags.FACTOR_LHS)
		ws.kkt_solver = ws.settings.kkt_solver_type(
							SparseMatrixCSC(ws.p.P),
							SparseMatrixCSC(ws.p.A),
							ws.settings.sigma,
							ws.ρvec)
	end
end
