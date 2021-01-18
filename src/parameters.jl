# set initial values of rhoVec
function set_rho_vec!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}
	m = ws.p.model_size[1]

	ws.ρ = ws.settings.rho
	ws.ρvec = ws.ρ * ones(T, m)

	apply_constraint_rho_scaling!(ws.ρvec, ws.row_ranges, ws.p.C.sets, ws.settings.RHO_MIN, ws.settings.RHO_EQ_OVER_RHO_INEQ)

	push!(ws.rho_updates, ws.ρ)
	return nothing
end


"Scale active constraints by `RHO_EQ_OVER_RHO_INEQ` and set loose constraints to `RHO_MIN`."
function apply_constraint_rho_scaling!(ρvec::AbstractVector{T}, row_ranges::Array{UnitRange{Int}, 1}, sets::Vector{AbstractConvexSet}, RHO_MIN::T, RHO_EQ_OVER_RHO_INEQ::T) where {T <: AbstractFloat}
	for (k, C) in enumerate(sets)
		if C isa ZeroSet
			for row in row_ranges[k]
				ρvec[row] *= RHO_EQ_OVER_RHO_INEQ
			end
		elseif C isa Nonnegatives
			row_offset = row_ranges[k].start - 1
			for (j, val) in enumerate(C.constr_type)
				if val == true
					ρvec[row_offset + j] = RHO_MIN
				end
			end

		elseif C isa Box
			row_offset = row_ranges[k].start - 1
			_rho_scale_box_bounds!(ρvec, C.constr_type, row_offset, RHO_MIN::T, RHO_EQ_OVER_RHO_INEQ::T)
		end
	end
	return nothing
end

function _rho_scale_box_bounds!(ρvec::AbstractVector{T}, constr_type::Vector{Int}, row_offset::Int, RHO_MIN::T, RHO_EQ_OVER_RHO_INEQ::T) where {T <: AbstractFloat}
	@inbounds for j = eachindex(constr_type)
		# loose inequality constraint
		if constr_type[j] == -1
			ρvec[j + row_offset] = RHO_MIN
		# equality constraint
		elseif constr_type[j] == 1
			ρvec[j + row_offset] *= RHO_EQ_OVER_RHO_INEQ
		end
	end
end


# adapt rhoVec based on residual ratio
function adapt_rho_vec!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}
	settings = ws.settings
	# compute normalized residuals based on the working variables (dont unscale)
	ignore_scaling = true
	r_prim::T, r_dual::T = calculate_residuals!(ws, ignore_scaling)
	max_norm_prim::T, max_norm_dual::T  = max_res_component_norm(ws, ignore_scaling)

	r_prim = r_prim / (max_norm_prim + T(1e-10))
	r_dual = r_dual / (max_norm_dual + T(1e-10))

	new_rho = ws.ρ * sqrt(r_prim / (r_dual + T(1e-10)))
	new_rho = min(max(new_rho, settings.RHO_MIN), settings.RHO_MAX)
	# only update rho if significantly different than current rho
	if (new_rho > settings.adaptive_rho_tolerance * ws.ρ) || (new_rho < (one(T) ./ settings.adaptive_rho_tolerance) * ws.ρ)
		update_rho_vec!(new_rho, ws)
	end
	return nothing
end

function update_rho_vec!(new_rho::T, ws::COSMO.Workspace{T}) where {T <: AbstractFloat}

	ws.ρ     = new_rho
	ws.ρvec .= new_rho

	apply_constraint_rho_scaling!(ws.ρvec, ws.row_ranges, ws.p.C.sets, ws.settings.RHO_MIN, ws.settings.RHO_EQ_OVER_RHO_INEQ)

	# log rho updates
	push!(ws.rho_updates, new_rho)

	if ws.settings.verbose_timing
		ws.times.factor_update_time += @elapsed update_rho!(ws.kkt_solver, ws.ρvec)
	else
		update_rho!(ws.kkt_solver, ws.ρvec)
	end

	return nothing
end
