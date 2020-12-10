#  Compute r = A * x + s - b in-place
function primal_kkt_condition!(r::AbstractVector{T}, A::AbstractMatrix{T}, x::AbstractVector{T}, s::SplitVector{T}, b::AbstractVector{T}) where {T <: AbstractFloat}
	# r = A * x
	mul!(r, A, x)
	@. r += s
	@. r -= b
	return nothing
end

# Compute r = P * x + q - A' * μ in-place
function dual_kkt_condition!(r::AbstractVector{T}, r_temp::AbstractVector{T}, P::AbstractMatrix{T}, x::AbstractVector{T}, q::AbstractVector{T}, A::AbstractMatrix{T}, μ::AbstractVector{T}) where {T <: AbstractFloat}
	mul!(r, P, x)
	@. r += q

	mul!(r_temp, A', μ)
	@. r -= r_temp
	return nothing
end

# Einv * x
unscale_primal!(x::AbstractVector{T}, Einv::AbstractMatrix{T}) where {T} = lmul!(Einv, x)

# cinv * Dinv * x
function unscale_dual!(x::AbstractVector{T}, Dinv::AbstractMatrix{T}, cinv) where {T}
		lmul!(Dinv, x)
		@. x *= cinv
		return nothing
end

function calculate_residuals!(ws::COSMO.Workspace,  IGNORESCALING_FLAG::Bool = false)
	# reuse vectors from temporary workspace
	r_prim = ws.utility_vars.vec_m
	r_dual = ws.utility_vars.vec_n
	r_temp = ws.utility_vars.vec_n2

	# r_prim = A * x + s - b
	primal_kkt_condition!(r_prim, ws.p.A, ws.vars.x, ws.vars.s, ws.p.b)

	# ∇f0 + ∑ νi ∇hi(x*) == 0 condition
	# r_dual = P * x + q - A' * μ
	dual_kkt_condition!(r_dual, r_temp, ws.p.P, ws.vars.x, ws.p.q, ws.p.A, ws.vars.μ)

	if (ws.settings.scaling != 0 && !IGNORESCALING_FLAG)
		# unscale primal residual: r_prim = Einv * ( A * x + s - b)
		unscale_primal!(r_prim, ws.sm.Einv)

		# unscale dual residual: r_dual = cinv * Dinv * (P * x + q - A' * μ)
		unscale_dual!(r_dual, ws.sm.Dinv, ws.sm.cinv)
	end

	return norm(r_prim, Inf), norm(r_dual, Inf)

end


function max_res_component_norm(ws::COSMO.Workspace, IGNORESCALING_FLAG::Bool = false)
		r_prim = ws.utility_vars.vec_m
		r_dual = ws.utility_vars.vec_n

		# Determine if the components inside max{ } have to be unscaled
		unscale = (ws.settings.scaling != 0 && !IGNORESCALING_FLAG)

		# Calculate: max { ||Einv * Ax ||_Inf, ||Einv * s||_Inf, ||Einv * b||_Inf }
		# Einv * A * x
		mul!(r_prim, ws.p.A, ws.vars.x)
		unscale && unscale_primal!(r_prim, ws.sm.Einv)
		max_norm_prim = norm(r_prim, Inf)

		# Einv * s
		@. r_prim = ws.vars.s
		unscale && unscale_primal!(r_prim, ws.sm.Einv)
		max_norm_prim = max.(max_norm_prim, norm(r_prim, Inf))

		# Einv * b
		@. r_prim = ws.p.b
		unscale && unscale_primal!(r_prim, ws.sm.Einv)
		max_norm_prim = max.(max_norm_prim, norm(r_prim, Inf))

		# Calculate: max { ||cinv * Dinv * P * x ||_Inf, ||cinv * Dinv * q||_Inf, ||cinv * Dinv * A' * μ||_Inf }
		# ||cinv * Dinv * P * x ||_Inf
		mul!(r_dual, ws.p.P, ws.vars.x)
		unscale && unscale_dual!(r_dual, ws.sm.Dinv, ws.sm.cinv)
		max_norm_dual = norm(r_dual, Inf)

		# ||cinv * Dinv * q||_Inf
		@. r_dual = ws.p.q
		unscale && unscale_dual!(r_dual, ws.sm.Dinv, ws.sm.cinv)
		max_norm_dual = max.(max_norm_dual, norm(r_dual, Inf))

		# ||cinv * Dinv * A' * μ||_Inf
		mul!(r_dual, ws.p.A', ws.vars.μ)
		unscale && unscale_dual!(r_dual, ws.sm.Dinv, ws.sm.cinv)
		max_norm_dual = max.(max_norm_dual, norm(r_dual, Inf))

	return max_norm_prim, max_norm_dual
end

function has_converged(ws::COSMO.Workspace{T}, r_prim::T, r_dual::T) where {T <: AbstractFloat}
	max_norm_prim, max_norm_dual = max_res_component_norm(ws)
	settings = ws.settings
	ϵ_prim = settings.eps_abs + settings.eps_rel * max_norm_prim
	ϵ_dual = settings.eps_abs + settings.eps_rel * max_norm_dual

	# check activation of accelerator
	COSMO.check_activation!(ws.accelerator, r_prim, r_dual, max_norm_prim, max_norm_dual)


	# if an optimal objective value was specified for the problem check if current solution is within specified accuracy
	obj_true_FLAG = true
	if !isnan(settings.obj_true)
		current_cost = calculate_cost!(ws.utility_vars.n, ws.vars.x, ws.p.P, ws.p.q, ws.sm.cinv[])
		obj_true_FLAG = abs(settings.obj_true - current_cost) <= settings.obj_true_tol
	end

	return (r_prim < ϵ_prim  && r_dual < ϵ_dual && obj_true_FLAG)
end

# cost = cinv *( 1/2 x' * P * x + q' x)
function calculate_cost!(temp::AbstractVector{T}, x::AbstractVector{T}, P::SparseMatrixCSC{T, Int}, q::AbstractVector{T}, cinv::T = one(T)) where {T <: AbstractFloat}
	# P x
	mul!(temp, P, x)
	return cinv * (T(0.5) * dot(temp, x) + dot(q, x))
end


