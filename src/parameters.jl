# set initial values of rhoVec
function set_rho_vec!(ws::COSMO.Workspace)
	m = ws.p.model_size[1]
	# nEQ = p.K.f
	# nINEQ = p.m - p.K.f
	ws.ρ = ws.settings.rho
	ws.ρvec = ws.ρ * ones(m)

	# scale ρ values that belong to equality constraints with a factor of 1e3
	set_ind = findall(x -> typeof(x) == COSMO.ZeroSet{Float64}, ws.p.C.sets)
	if length(set_ind) > 0
		row_ind = COSMO.get_set_indices(ws.p.C.sets)
		for (i, rows) in enumerate(row_ind[set_ind])
			ws.ρvec[rows] *= 1e3
		end
	end
	push!(ws.Info.rho_updates, ws.ρ)
	return nothing
end


# adapt rhoVec based on residual ratio
function adapt_rho_vec!(ws::COSMO.Workspace)
	settings = ws.settings
	# compute normalized residuals based on the working variables (dont unscale)
	ignore_scaling = true
	r_prim::Float64, r_dual::Float64 = calculate_residuals(ws, ignore_scaling)
	max_norm_prim::Float64, max_norm_dual::Float64  = max_res_component_norm(ws, ignore_scaling)

	r_prim = r_prim / (max_norm_prim + 1e-10)
	r_dual = r_dual / (max_norm_dual + 1e-10)

	new_rho = ws.ρ * sqrt(r_prim / (r_dual + 1e-10))
	new_rho = min(max(new_rho, settings.RHO_MIN), settings.RHO_MAX)
	# only update rho if significantly different than current rho
	if (new_rho > settings.adaptive_rho_tolerance * ws.ρ) || (new_rho < (1 ./ settings.adaptive_rho_tolerance) * ws.ρ)
		update_rho_vec!(new_rho, ws)
	end
	return nothing
end

function update_rho_vec!(new_rho::Float64, ws::COSMO.Workspace)

	ws.ρ     = new_rho
	ws.ρvec .= new_rho

	# scale ρ values that belong to equality constraints with a factor of 1e3
	set_ind = findall(x -> typeof(x) == COSMO.ZeroSet{Float64}, ws.p.C.sets)
	if length(set_ind) > 0
		row_ind = COSMO.get_set_indices(ws.p.C.sets)
		for (i, rows) in enumerate(row_ind[set_ind])
			ws.ρvec[rows] *= 1e3
		end
	end

	# log rho updates to info variable
	push!(ws.Info.rho_updates, new_rho)

	if ws.settings.verbose_timing
		ws.times.factor_time += @elapsed update_rho!(ws.kkt_solver,ws.ρvec)
	else
		update_rho!(ws.kkt_solver,ws.ρvec)
	end

	return nothing
end
