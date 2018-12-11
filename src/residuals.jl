function calculate_residuals(ws::COSMO.Workspace,  IGNORESCALING_FLAG::Bool = false)

	if (ws.settings.scaling != 0 && !IGNORESCALING_FLAG)
		r_prim = norm(ws.sm.Einv * (ws.p.A * ws.vars.x + ws.vars.s - ws.p.b), Inf)
		# ∇f0 + ∑ νi ∇hi(x*) == 0 condition
		r_dual = norm(ws.sm.cinv[] * ws.sm.Dinv * (ws.p.P * ws.vars.x + ws.p.q - ws.p.A' * ws.vars.μ), Inf)
	end
	if (ws.settings.scaling == 0 || IGNORESCALING_FLAG )
		r_prim = norm(ws.p.A * ws.vars.x + ws.vars.s - ws.p.b, Inf)
		r_dual = norm(ws.p.P * ws.vars.x + ws.p.q - ws.p.A' * ws.vars.μ, Inf)
	end
	# FIXME: Why is it -A'μ ?
	return r_prim,r_dual
end

function max_res_component_norm(ws::COSMO.Workspace, IGNORESCALING_FLAG::Bool = false)
	if (ws.settings.scaling != 0 && !IGNORESCALING_FLAG)
		max_norm_prim = max.(norm(ws.sm.Einv * (ws.p.A * ws.vars.x), Inf), norm(ws.sm.Einv * ws.vars.s, Inf), norm(ws.sm.Einv * ws.p.b, Inf))
		max_norm_dual = max.(norm(ws.sm.cinv .* (ws.sm.Dinv * (ws.p.P * ws.vars.x)), Inf), norm(ws.sm.cinv .* (ws.sm.Dinv * ws.p.q), Inf), norm(ws.sm.cinv .* (ws.sm.Dinv * (ws.p.A' * ws.vars.μ)), Inf))
	end
	if (ws.settings.scaling == 0 || IGNORESCALING_FLAG)
		max_norm_prim = max.(norm(ws.p.A * ws.vars.x, Inf), norm(ws.vars.s, Inf), norm(ws.p.b, Inf))
		max_norm_dual = max.(norm(ws.p.P * ws.vars.x, Inf), norm(ws.p.q, Inf), norm(ws.p.A' * ws.vars.μ, Inf))
	end
	return max_norm_prim, max_norm_dual
end

function has_converged(ws::COSMO.Workspace, r_prim::Float64, r_dual::Float64)
	max_norm_prim, max_norm_dual = max_res_component_norm(ws)
	settings = ws.settings
	ϵ_prim = settings.eps_abs + settings.eps_rel * max_norm_prim
	ϵ_dual = settings.eps_abs + settings.eps_rel * max_norm_dual

	# if an optimal objective value was specified for the problem check if current solution is within specified accuracy
	obj_true_FLAG = true
	if !isnan(settings.obj_true)
		current_cost = ws.sm.cinv[] * (1/2 * ws.vars.x' * ws.p.P * ws.vars.x + ws.p.q' * ws.vars.x)[1]
		obj_true_FLAG = abs(settings.obj_true - current_cost) <= settings.obj_true_tol
	end

	return (r_prim < ϵ_prim  && r_dual < ϵ_dual && obj_true_FLAG)
end
