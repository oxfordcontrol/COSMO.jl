# sup_{z in K_tilde_b = {-K} x {b} } <z,δy> = { <y,b> ,if y in Ktilde_polar
#                                                 +∞   ,else}

function support_function(y, ws)
	if in_dual(y, ws.p.C, ws.settings.eps_prim_inf)
		return dot(y, ws.p.b)
	else
		return Inf
	end

end

function is_primal_infeasible(δy, ws)
	settings = ws.settings
	# calculate unscaled norm of δy
	norm_δy = scaled_norm(ws.sm.E, δy, Inf)::eltype(δy)

	# make sure norm is unequal to zero before continuing
	if norm_δy > settings.eps_prim_inf

		# test condition A'δy = 0
		A_δy = ws.p.A' * δy
		A_δy = ws.sm.Dinv * A_δy

		if norm(A_δy, Inf) / norm_δy <= settings.eps_prim_inf
			# test condition S_K(δy) < 0
			unit_δy_split = SplitVector(δy / norm_δy, ws.p.C)
			sF = support_function(unit_δy_split, ws)
			if sF <= settings.eps_prim_inf
				return true
			end
		end
	end
	return false
end


function is_dual_infeasible(δx, ws)
	settings = ws.settings
	# calculate unscaled norm of δx
	norm_δx = scaled_norm(ws.sm.D, δx, Inf)::eltype(δx)

	if norm_δx > settings.eps_dual_inf
		# test condition <q,δx> < 0
		if  dot(ws.p.q,δx) / (norm_δx * ws.sm.c[]) < -settings.eps_dual_inf
			# test condition Pδx == 0
			P_δx = ws.p.P * δx
			P_δx = ws.sm.Dinv * P_δx
			if norm(P_δx, Inf) / (norm_δx * ws.sm.c[]) <= settings.eps_dual_inf

				# test condition Ax in Ktilde_b∞
				A_δx = ws.p.A * δx
				A_δx = ws.sm.Einv * A_δx
				A_δx_split = SplitVector(A_δx / norm_δx, ws.p.C)
				if in_recc(A_δx_split, ws.p.C, settings.eps_dual_inf)
					return true
				end
			end
		end
	end
	return false
end
