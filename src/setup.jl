function setup!(ws::COSMO.Workspace)
	# scale problem data
	if ws.settings.scaling != 0
		scale_ruiz!(ws)
	end

	set_rho_vec!(ws)

	# create a KKT solver object populated with our data
	if(ws.flags.FACTOR_LHS)
		if ws.settings.verbose_timing
			ws.times.factor_time += @elapsed _make_kkt_solver(ws)
		else
			_make_kkt_solver(ws)
		end

	end
end


function _make_kkt_solver(ws)

	ws.kkt_solver = ws.settings.kkt_solver_type(
						SparseMatrixCSC(ws.p.P),
						SparseMatrixCSC(ws.p.A),
						ws.settings.sigma,
						ws.ρvec)

end
