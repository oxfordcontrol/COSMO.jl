function _make_kkt_solver!(ws::COSMO.Workspace)

	ws.kkt_solver = ws.settings.kkt_solver(SparseMatrixCSC(ws.p.P),
						SparseMatrixCSC(ws.p.A),
						ws.settings.sigma,
						ws.ρvec)

end

function setup!(ws::COSMO.Workspace)

  allocate_set_memory!(ws)
	# scale problem data
	if ws.settings.scaling != 0
		scale_ruiz!(ws)
	end

	set_rho_vec!(ws)

	# create a KKT solver object populated with our data
	if(ws.flags.FACTOR_LHS)
		if ws.settings.verbose_timing
			ws.times.factor_time += @elapsed _make_kkt_solver!(ws)
		else
			_make_kkt_solver!(ws)
		end

	end
end

function allocate_set_memory!(ws)
  for (k, C) in enumerate(ws.p.C.sets)
    allocate_memory!(C)
  end
end
