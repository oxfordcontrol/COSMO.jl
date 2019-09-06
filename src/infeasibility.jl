function is_primal_infeasible!(δy, ws)

	settings = ws.settings
	# calculate unscaled norm of δy
	norm_δy = scaled_norm(ws.sm.E, δy.data, Inf)::eltype(δy)

	# make sure norm is unequal to zero before continuing
	if norm_δy > settings.eps_prim_inf
		# test condition A'δy = 0
		#A_δy = ws.p.A' * δy
		A_δy  = ws.utility_vars.vec_n
		mul!(A_δy, ws.p.A', δy.data)
		# unscale: A_δy = Dinv * A_δy
		lmul!(ws.sm.Dinv, A_δy)


		if norm(A_δy, Inf) <= settings.eps_prim_inf * norm_δy
			# scale δy to unit size
			@. δy *= (-1 / norm_δy)
			δyt_b = dot(δy, ws.p.b)
			# notice the in-place function. This is faster as we are using minus_unit_δy as workspace
			sF = support_function!(δy, ws.p.C, settings.eps_prim_inf) - δyt_b
			if sF <= settings.eps_prim_inf
				return true
			end
		end
	end
	return false
end


function is_dual_infeasible!(δx, ws)
	settings = ws.settings
	# calculate unscaled norm of δx
	norm_δx = scaled_norm(ws.sm.D, δx, Inf)::eltype(δx)

	if norm_δx > settings.eps_dual_inf
		# test condition <q,δx> < 0
		if  dot(ws.p.q, δx) / (norm_δx * ws.sm.c[]) < -settings.eps_dual_inf

			# test condition Pδx == 0
			P_δx = ws.utility_vars.vec_n
			# P_δx = P * δx
			mul!(P_δx, ws.p.P, δx)

			# unscale P_δx = Dinv * P_δx
			lmul!(ws.sm.Dinv, P_δx)

			if norm(P_δx, Inf) / (norm_δx * ws.sm.c[]) <= settings.eps_dual_inf

				# test condition Ax in Ktilde_b∞
				A_δx = ws.utility_vars.vec_m
				mul!(A_δx, ws.p.A, δx)

				# unscale A_δx = Einv * P_δx
				lmul!(ws.sm.Einv, A_δx)

				#normalize A_δx and create a SplitVector
				A_δx .*= 1/norm_δx
				A_δx_split = SplitVector(A_δx, ws.p.C)
				if in_pol_recc!(A_δx_split, ws.p.C, settings.eps_dual_inf)
					return true
				end
			end
		end
	end
	return false
end
