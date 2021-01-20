function _make_kkt_solver!(ws::COSMO.Workspace)

	ws.kkt_solver = ws.settings.kkt_solver(SparseMatrixCSC(ws.p.P),
						SparseMatrixCSC(ws.p.A),
						ws.settings.sigma,
						ws.ρvec)
end


function _make_accelerator!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}
  # if the user passed in a custom AbstractAccelerator, e.g. with different value for memory, don't change it
	m, n = ws.p.model_size
	ws.accelerator = ws.settings.accelerator(m + n)
	ws.accelerator_safeguarding = ws.settings.safeguarding
	ws.safeguarding_tol = ws.settings.safeguarding_tol
end



function setup!(ws::COSMO.Workspace)

  	allocate_set_memory!(ws)
	# scale problem data
	if ws.settings.scaling != 0 && !ws.states.IS_SCALED
		if ws.settings.verbose_timing
			ws.times.scaling_time = @elapsed scale_ruiz!(ws)
		else
			scale_ruiz!(ws)
		end
		ws.states.IS_SCALED = true
	else
		# even if the problem data is already scaled, we still have to scale our (potentially warm-started) variables
		scale_variables!(ws.vars.x, ws.vars.μ, ws.vars.s, ws.sm.Dinv, ws.sm.Einv, ws.sm.E, ws.sm.c)
		ws.times.scaling_time = 0.
	end

	# construct a set index --> row range map
	ws.row_ranges = get_set_indices(ws.p.C.sets)
	classify_constraints!(ws)

	# only set ρ to settings value if this is the first time solving, otherwise keep ρ in sync with value used in KKT-factorisation
	if !ws.states.IS_OPTIMIZED
		set_rho_vec!(ws)
	end

	# instantiate accelerator
	if !ws.states.KKT_FACTORED
		_make_accelerator!(ws)
	else
		CA.restart!(ws.accelerator)
		ws.accelerator.activated = false
	end


	# create a KKT solver object populated with our data
	if !ws.states.KKT_FACTORED
		if ws.settings.verbose_timing
			ws.times.init_factor_time = @elapsed _make_kkt_solver!(ws)
		else
			_make_kkt_solver!(ws)
		end
		ws.states.KKT_FACTORED = true
	end


end

function allocate_set_memory!(ws::COSMO.Workspace)
	if !ws.states.IS_OPTIMIZED
		for (k, C) in enumerate(ws.p.C.sets)
			allocate_memory!(C)
		end
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
