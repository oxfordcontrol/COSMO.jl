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
		if ws.settings.verbose_timing
			ws.times.scaling_time = @elapsed scale_ruiz!(ws)
		else
			scale_ruiz!(ws)
		end
	else
		ws.times.scaling_time = 0.
	end

	# construct a set index --> row range map
	ws.row_ranges = get_set_indices(ws.p.C.sets)
	classify_constraints!(ws)

	set_rho_vec!(ws)

	# create a KKT solver object populated with our data
	if(ws.flags.FACTOR_LHS)
		if ws.settings.verbose_timing
			ws.times.init_factor_time = @elapsed _make_kkt_solver!(ws)
		else
			_make_kkt_solver!(ws)
		end

	end
end

function allocate_set_memory!(ws::COSMO.Workspace)
  for (k, C) in enumerate(ws.p.C.sets)
    allocate_memory!(C)
  end
end

"For inequality and box constraints, check if any constraints are obviously loose or tight."
function classify_constraints!(ws::COSMO.Workspace)
	for (k, C) in enumerate(ws.p.C.sets)
		if C isa Nonnegatives
			b = view(ws.p.b, ws.row_ranges[k])
			classify_constraints!(C.constr_type, b, ws.settings.COSMO_INFTY, ws.settings.MIN_SCALING)
		elseif C isa Box
			classify_box_constraints!(C, ws.settings.COSMO_INFTY, ws.settings.MIN_SCALING, ws.settings.RHO_TOL)
		end
	end
	return nothing
end
